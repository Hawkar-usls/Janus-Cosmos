#!/usr/bin/env python3
"""R2E: prelocked LOVE high-time-resolution Parkes raw file qualification.
Only official DAP TAP metadata, never raw spectral power and never signal claim.
"""
import hashlib,json,math,pathlib
import osiris_spear_r2d_casda_psrda as r2d
P=pathlib.Path("research/seti/osiris_spear_r2e_love_protocol.json")
S=json.loads(P.read_text())
OUT=pathlib.Path("runtime/osiris_spear_r2e");OUT.mkdir(parents=True,exist_ok=True)
SOURCE=r2d.ROUTES["PSRDA_DAP"]
def save(f,x):(OUT/f).write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
def num(x):
    try:
        z=float(x)
        return z if math.isfinite(z) else None
    except (ValueError,TypeError):return None
def ang(a,b,c,d):return r2d.old.angle(a,b,c,d)
def main():
    meta=r2d.schema("PSRDA_DAP")
    assert meta["schema_state"]=="OK"
    names={x["name"] for x in meta["columns"]}
    fields=[x for x in S["inspect_metadata_fields"] if x in names]
    assert {"projid","filename","ra_angle","dec_angle","obsfreq","obs_mode","obstype"}.issubset(fields)
    result=[]
    ra,dec=S["fixed_target"]["ra_deg"],S["fixed_target"]["dec_deg"]
    locked=S["preselected_original_files"]
    for item in locked:
        file=item["filename"].replace("'","''")
        q=f"SELECT TOP 20 {','.join(fields)} FROM psrda.pulsar_obs WHERE filename='{file}' AND projid='{item['project']}'"
        receipt,rows=r2d.tap(SOURCE,q,20)
        identities=[]
        for row in rows:
            ra0=num(row.get("ra_angle"));de0=num(row.get("dec_angle"))
            dist=ang(ra,dec,ra0,de0) if ra0 is not None and de0 is not None else None
            other={f:r2d.safe(row.get(f)) for f in fields}
            other.update({"distance_to_frozen_LOVE_deg":dist,
              "matches_prefrozen_selection":abs(dist-item["angular_offset_deg"])<.0001 if dist is not None else False,
              "archive_raw_SEARCH":str(row.get("obstype","")).lower()=="raw" and str(row.get("obs_mode","")).upper()=="SEARCH",
              "claimed_receiver_verified_by_TCS":False,
              "file_byte_hash_verified":False})
            identities.append(other)
        r={"file":item["filename"],"project":item["project"],"lookup":receipt,
           "metadata_rows":len(rows),"identities":identities}
        result.append(r)
        print("OSIRIS_R2E_RAW_IDENTITY",json.dumps({"file":item["filename"],"query_ok":receipt["state"],
          "rows":len(rows),"distance_deg":identities[0].get("distance_to_frozen_LOVE_deg") if identities else None,
          "frontend":identities[0].get("frontend") if identities else None,"backend":identities[0].get("backend") if identities else None,
          "filesize":identities[0].get("filesize") if identities else None,
          "collection_fedora_pid":identities[0].get("collection_fedora_pid") if identities else None},ensure_ascii=False),flush=True)
    save("file_identity_and_metadata.json",{"data":result,"original_bytes_verified":0,"not_raw_signal_analysis":True})
    groups=[]
    # Prespecified same-session *metadata* panel; NOT true OFF observations.
    for proj,prefix in (("P892","bpsr170815_%"),("P360","BJ0093_%")):
        q=f"SELECT TOP 200 {','.join(fields)} FROM psrda.pulsar_obs WHERE projid='{proj}' AND filename LIKE '{prefix}'"
        receipt,rows=r2d.tap(SOURCE,q,200)
        summaries=[]
        for row in rows:
            ra0=num(row.get("ra_angle"));de0=num(row.get("dec_angle"))
            dist=ang(ra,dec,ra0,de0) if ra0 is not None and de0 is not None else None
            summaries.append({"filename":row.get("filename"),"date_obs":row.get("date_obs"),
                 "obsfreq":row.get("obsfreq"),"obstype":row.get("obstype"),"obs_mode":row.get("obs_mode"),
                 "angular_distance_LOVE_deg":dist,
                 "frontend":row.get("frontend"),"backend":row.get("backend")})
        groups.append({"project":proj,"prefix":prefix,"query_receipt":receipt,"rows_returned":len(rows),
           "truncated":receipt.get("maybe_truncated",False),"examples":summaries[:120],
           "matching_OFF_present_and_qualified":False})
        print("OSIRIS_R2E_SAME_EPOCH",json.dumps({"proj":proj,"state":receipt["state"],"rows":len(rows),
          "truncated":receipt.get("maybe_truncated",False)},ensure_ascii=False),flush=True)
    save("epoch_control_survey.json",{"groups":groups,
         "multibeam_or_temporal_nearby_not_independent_OFF_without_actual_beam":True})
    good=[]
    for item in result:
        for x in item["identities"]:
            if x["matches_prefrozen_selection"] and x["archive_raw_SEARCH"]:good.append(x)
    report={"stage":"OSIRIS_SPEAR_R2E","status":
      "ORIGINAL_RAW_SEARCH_METADATA_IDENTITIES_FOUND__HOLD_NATIVE_BYTES_POINTING_AND_OFF" if good else
      "RAW_SEARCH_SOURCE_METADATA_UNQUALIFIED__HOLD",
      "fixed_target":"LOVE_G1_007","raw_search_candidates_preselected":len(locked),
      "verified_raw_search_identity_rows":len(good),
      "files":[{"filename":x["filename"],"project":x["projid"],
        "angular_sep_deg":x["distance_to_frozen_LOVE_deg"],"front_end":x.get("frontend"),
        "back_end":x.get("backend"),"file_size":x.get("filesize"),
        "collection_id":x.get("collection_fedora_pid"),"obsbw":x.get("obsbw"),
        "obsnchan":x.get("obsnchan"),"timestamp":x.get("date_obs")} for x in good],
      "native_raw_files_downloaded":0,"source_pointing_antenna_TCS_verified":False,
      "beam_actual_receiver_certified":False,"matched_ONOFF_qualified":False,
      "scientific_frequency_axis_verified":False,"actual_JANUS_detector_run":False,
      "false_technosignature_claim":False,"non_detection_claim":False,
      "no_production_model_change":True,"R42B_external_gate":"HOLD",
      "next":"Confirm PSRFITS header per-subintegration actual telescope pointing, beam index and gain, original native files and SHA256, and same-session physically independent controls. A nominal 14-14.5 arcmin FWHM at 1.4GHz makes 6.06/6.20 arcmin near the nominal 7-7.25 arcmin half-power radius only conditionally."}
    save("r2e_gate_report.json",report)
    print("OSIRIS_R2E_RESULT",json.dumps(report,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
