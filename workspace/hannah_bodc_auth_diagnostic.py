#!/usr/bin/env python3
"""Non-secret transport diagnostic for Hannah/BODC archive access.

No requester credentials are used. The diagnostic records public SSH auth
methods and tests conventional anonymous FTP using an invalid generic contact
address. It never attempts to infer or brute-force credentials.
"""
from __future__ import annotations

import ftplib
import hashlib
import json
import socket
from datetime import datetime, timezone
from pathlib import Path

import paramiko

HOST = "livftp.noc.ac.uk"
GENERIC_ANON = "janus-probe@example.invalid"
REQUESTS = ("BODCREQ-9408", "BODCREQ-9406")


def ssh_diag() -> dict:
    out = {"port": 22, "reachable": False, "allowed_auth_methods": []}
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
    finally:
        t.close()
    return out


def safe_nlst(ftp: ftplib.FTP, path: str | None = None) -> list[str]:
    try:
        return sorted(ftp.nlst(path)) if path else sorted(ftp.nlst())
    except Exception:
        return []


def ftp_plain_diag() -> dict:
    out = {"port": 21, "tls": False, "reachable": False, "anonymous_generic_login": False}
    ftp = ftplib.FTP(timeout=20)
    try:
        out["banner"] = ftp.connect(HOST, 21, timeout=20)
        out["reachable"] = True
        out["login_reply"] = ftp.login("anonymous", GENERIC_ANON)
        out["anonymous_generic_login"] = True
        out["pwd_after_login"] = ftp.pwd()
        root_names = safe_nlst(ftp)
        out["root_listing_names"] = root_names[:200]
        out["root_listing_truncated"] = len(root_names) > 200

        candidate_templates = [
            "/{req}",
            "{req}",
            "/data/{req}",
            "data/{req}",
            "/bodc/data/{req}",
            "bodc/data/{req}",
            "/bodc/bodc/data/{req}",
            "bodc/bodc/data/{req}",
        ]
        located = {}
        for req in REQUESTS:
            trials = []
            found = None
            for template in candidate_templates:
                candidate = template.format(req=req)
                try:
                    original = ftp.pwd()
                    ftp.cwd(candidate)
                    found = ftp.pwd()
                    names = safe_nlst(ftp)
                    trials.append({"candidate": candidate, "ok": True, "resolved_pwd": found})
                    located[req] = {
                        "resolved_path": found,
                        "listing_names": names[:200],
                        "listing_truncated": len(names) > 200,
                    }
                    ftp.cwd(original)
                    break
                except Exception as exc:
                    trials.append({"candidate": candidate, "ok": False, "failure": str(exc)})
                    try:
                        ftp.cwd(out["pwd_after_login"])
                    except Exception:
                        pass
            if found is None:
                located[req] = {"resolved_path": None, "trials": trials}
        out["requests"] = located
    except Exception as exc:
        out["failure_type"] = type(exc).__name__
        out["failure"] = str(exc)
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
        "schema": "janus.cosmos.cousteau.hannah_bodc.auth_transport_diagnostic.v2",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "uses_requester_credentials": False,
        "bruteforce": False,
        "host": HOST,
    }
    try:
        result["ssh_sftp"] = ssh_diag()
    except Exception as exc:
        result["ssh_sftp"] = {"reachable": False, "failure_type": type(exc).__name__, "failure": str(exc)}
    result["ftp_plain"] = ftp_plain_diag()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
