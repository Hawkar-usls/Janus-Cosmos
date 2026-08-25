#!/usr/bin/env python3
"""Non-secret transport diagnostic for Hannah/BODC archive access.

No requester credentials are used. The diagnostic records public SSH auth
methods and tests conventional anonymous FTP/FTPS using an invalid generic
contact address. It never attempts to infer or brute-force credentials.
"""
from __future__ import annotations

import ftplib
import hashlib
import json
import socket
import ssl
from datetime import datetime, timezone
from pathlib import Path

import paramiko

HOST = "livftp.noc.ac.uk"
GENERIC_ANON = "janus-probe@example.invalid"
REMOTE = "/bodc/bodc/data/BODCREQ-9408"


def ssh_diag() -> dict:
    out = {"port": 22, "reachable": False, "allowed_auth_methods": [], "interactive_prompts": []}
    sock = socket.create_connection((HOST, 22), timeout=15)
    t = paramiko.Transport(sock)
    try:
        t.start_client(timeout=15)
        out["reachable"] = True
        key = t.get_remote_server_key()
        out["host_key"] = {
            "type": key.get_name(),
            "sha256_hex": hashlib.sha256(key.asbytes()).hexdigest(),
        }
        try:
            t.auth_none("anonymous")
            out["none_auth_succeeded"] = True
        except paramiko.BadAuthenticationType as exc:
            out["none_auth_succeeded"] = False
            out["allowed_auth_methods"] = sorted(exc.allowed_types)
        except paramiko.AuthenticationException:
            out["none_auth_succeeded"] = False

        if "keyboard-interactive" in out["allowed_auth_methods"]:
            seen = []
            def handler(title, instructions, prompts):
                seen.append({
                    "title": title,
                    "instructions": instructions,
                    "prompts": [{"prompt": p, "echo": bool(e)} for p, e in prompts],
                })
                return [""] * len(prompts)
            try:
                t.auth_interactive("anonymous", handler)
                out["empty_interactive_succeeded"] = True
            except Exception as exc:
                out["empty_interactive_succeeded"] = False
                out["interactive_failure_type"] = type(exc).__name__
            out["interactive_prompts"] = seen
    finally:
        t.close()
    return out


def ftp_diag(tls: bool) -> dict:
    out = {"port": 21, "tls": tls, "reachable": False, "anonymous_generic_login": False}
    cls = ftplib.FTP_TLS if tls else ftplib.FTP
    ftp = cls(timeout=20)
    try:
        banner = ftp.connect(HOST, 21, timeout=20)
        out["reachable"] = True
        out["banner"] = banner
        if tls:
            ftp.auth()
            ftp.prot_p()
        try:
            reply = ftp.login("anonymous", GENERIC_ANON)
            out["anonymous_generic_login"] = True
            out["login_reply"] = reply
            try:
                ftp.cwd(REMOTE)
                out["request_path_accessible"] = True
                names = ftp.nlst()
                out["request_listing_names"] = sorted(names)[:100]
                out["request_listing_truncated"] = len(names) > 100
            except Exception as exc:
                out["request_path_accessible"] = False
                out["path_failure_type"] = type(exc).__name__
                out["path_failure"] = str(exc)
        except Exception as exc:
            out["login_failure_type"] = type(exc).__name__
            out["login_failure"] = str(exc)
    except Exception as exc:
        out["connection_failure_type"] = type(exc).__name__
        out["connection_failure"] = str(exc)
    finally:
        try:
            ftp.quit()
        except Exception:
            try:
                ftp.close()
            except Exception:
                pass
    return out


def main() -> int:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()
    result = {
        "schema": "janus.cosmos.cousteau.hannah_bodc.auth_transport_diagnostic.v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "uses_requester_credentials": False,
        "bruteforce": False,
        "host": HOST,
    }
    try:
        result["ssh_sftp"] = ssh_diag()
    except Exception as exc:
        result["ssh_sftp"] = {"reachable": False, "failure_type": type(exc).__name__, "failure": str(exc)}
    result["ftp_plain"] = ftp_diag(False)
    result["ftp_tls"] = ftp_diag(True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
