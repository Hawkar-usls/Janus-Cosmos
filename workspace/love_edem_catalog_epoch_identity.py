#!/usr/bin/env python3
"""LIVE LOVE/EDEM catalog-epoch identity gate.

Resolves exact/reference epochs and catalog-specific positional uncertainty for
already-frozen source IDs, then asks whether those independent detections are
compatible with the corresponding Gaia DR3 stellar worldline.
"""
from __future__ import annotations
import argparse,csv,io,json,math
from datetime import datetime,timezone
from pathlib import Path
from typing import Any,Dict,Iterable
import requests
from janus_catalog_epoch_identity import CatalogMeasurement,conservative_unknown_corr_bound_mas2,covariance_from_cosigma_mas2,ellipse_covariance_mas2,evaluate_identity,jd_to_mjd
from love_edem_epoch_worldline_test import query_gaia_live,state_from_features

SCHEMA="janus.cosmos.love-edem.catalog-epoch-identity.v1"
SDSS_SQL_ENDPOINT="https://skyserver.sdss.org/dr16/SkyServerWS/SearchTools/SqlSearch"
PS1_SYSTEMATIC_FLOOR_MAS=20.0
GROUPS={
 "LOVE_GAIA_WISE_2MASS":{"target":"LOVE","gaia_source_id":"6163586620213012352","center":[204.29573215827,-36.77911588007],"catalogs":{"ALLWISE":"J133710.98-364644.9","2MASS_PSC":"13371100-3646446"}},
 "EDEM_GAIA_WISE_SDSS_PS1":{"target":"EDEM_SEARCH_CENTER_ZP","gaia_source_id":"699051350998534656","center":[139.21802643661,30.26276387854],"catalogs":{"ALLWISE":"J091652.31+301545.8","PANSTARRS_DR1":"144311392180445928","SDSS_DR16":"1237664320537362761"}},
}
CATALOG_META={
 "ALLWISE":{"vizier":"II/328/allwise","id_columns":["AllWISE"],"constraint_names":["AllWISE"]},
 "2MASS_PSC":{"vizier":"II/246/out","id_columns":["2MASS","_2MASS"],"constraint_names":["2MASS","_2MASS"]},
 "PANSTARRS_DR1":{"vizier":"II/349/ps1","id_columns":["objID"],"constraint_names":["objID"]},
 "SDSS_DR16":{"vizier":"V/154/sdss16","id_columns":["objID"],"constraint_names":["objID"]},
}

def _jsonable(v:Any)->Any:
 if v is None:return None
 try:
  if getattr(v,"mask",False):return None
 except Exception:pass
 if isinstance(v,(str,int,float,bool)):return v
 try:return v.item()
 except Exception:return str(v)

def _get(d:Dict[str,Any],*names:str)->Any:
 low={str(k).lower():k for k in d}
 for name in names:
  k=low.get(name.lower())
  if k is not None:
   v=d[k]
   if v is not None and str(v).strip() not in ("","--","nan"):return v
 return None

def _float(d:Dict[str,Any],*names:str)->float|None:
 v=_get(d,*names)
 if v is None:return None
 try:
  x=float(v);return x if math.isfinite(x) else None
 except Exception:return None

def _id_matches(row:Dict[str,Any],expected:str,id_columns:Iterable[str])->bool:
 e=str(expected).strip().replace("WISEA ","").replace("WISE ","")
 for c in id_columns:
  v=_get(row,c)
  if v is not None and str(v).strip().replace("WISEA ","").replace("WISE ","")==e:return True
 return False

def query_vizier_exact(catalog_key:str,source_id:str,center:list[float])->Dict[str,Any]:
 from astroquery.vizier import Vizier
 import astropy.units as u
 from astropy.coordinates import SkyCoord
 cfg=CATALOG_META[catalog_key];viz=Vizier(columns=["**"],row_limit=100);errors=[]
 for constraint in cfg["constraint_names"]:
  try:
   for table in viz.query_constraints(catalog=cfg["vizier"],**{constraint:str(source_id)}):
    for row in table:
     snap={n:_jsonable(row[n]) for n in table.colnames}
     if _id_matches(snap,source_id,cfg["id_columns"]):snap["_janus_query_mode"]=f"EXACT_CONSTRAINT:{constraint}";return snap
  except Exception as exc:errors.append(f"{constraint}:{type(exc).__name__}:{exc}")
 try:
  c=SkyCoord(center[0]*u.deg,center[1]*u.deg,frame="icrs")
  for table in viz.query_region(c,radius=2*u.arcsec,catalog=cfg["vizier"]):
   for row in table:
    snap={n:_jsonable(row[n]) for n in table.colnames}
    if _id_matches(snap,source_id,cfg["id_columns"]):snap["_janus_query_mode"]="POSITION_RECOVERY_WITH_EXACT_ID_ACCEPTANCE";snap["_janus_constraint_errors"]=errors;return snap
 except Exception as exc:errors.append(f"region:{type(exc).__name__}:{exc}")
 raise RuntimeError(f"EXACT_SOURCE_NOT_RETURNED:{catalog_key}:{source_id}:{errors}")

def allwise_measurement(source_id:str,row:Dict[str,Any])->CatalogMeasurement:
 ra,dec=_float(row,"RA_pm"),_float(row,"DE_pm");era,edec,co=_float(row,"e_RA_pm","sigra_pm"),_float(row,"e_DE_pm","sigdec_pm"),_float(row,"cosig_pm","sigradec_pm")
 if None in (ra,dec,era,edec,co):raise ValueError("ALLWISE_PM_REFERENCE_SOLUTION_INCOMPLETE")
 cov=covariance_from_cosigma_mas2(era,edec,co)
 return CatalogMeasurement("ALLWISE",source_id,ra,dec,55400.0,cov.tolist(),"MEASURED_FULL_2D_COVARIANCE_FROM_ALLWISE_PM_COSIGMA","EXACT_CATALOG_REFERENCE_EPOCH_MJD_55400","ALLWISE_MOTION_FIT_POSITION_AT_STANDARD_REFERENCE_EPOCH",{"catalog_id":"II/328/allwise","query_mode":row.get("_janus_query_mode"),"raw_epoch_evidence":{"RA_pm":ra,"DE_pm":dec,"reference_mjd":55400.0},"raw_uncertainty":{"e_RA_pm_arcsec":era,"e_DE_pm_arcsec":edec,"cosig_pm_arcsec":co},"stationary_position_not_used":True})

def twomass_measurement(source_id:str,row:Dict[str,Any])->CatalogMeasurement:
 ra,dec=_float(row,"RAJ2000"),_float(row,"DEJ2000");maj,minor,pa=_float(row,"errMaj"),_float(row,"errMin"),_float(row,"errPA");raw_jd=_float(row,"JD","jdate")
 if None in (ra,dec,maj,minor,pa,raw_jd):raise ValueError("2MASS_POSITION_EPOCH_OR_ELLIPSE_INCOMPLETE")
 mjd=jd_to_mjd(raw_jd);cov=ellipse_covariance_mas2(maj,minor,pa)
 return CatalogMeasurement("2MASS_PSC",source_id,ra,dec,mjd,cov.tolist(),"MEASURED_FULL_2D_ERROR_ELLIPSE","EXACT_OBSERVATION_EPOCH_FROM_2MASS_JD","2MASS_PSC_POSITION_AT_SCAN_OBSERVATION_EPOCH",{"catalog_id":"II/246/out","query_mode":row.get("_janus_query_mode"),"raw_epoch_field":"JD","raw_epoch_value":raw_jd,"interpreted_mjd":mjd,"date_precision_contract":"2MASS_JDATE_DOCUMENTED_TO_APPROX_PLUS_MINUS_30_SECONDS","raw_error_ellipse":{"major_arcsec":maj,"minor_arcsec":minor,"pa_deg_east_of_north":pa}})

def ps1_measurement(source_id:str,row:Dict[str,Any])->CatalogMeasurement:
 ra,dec=_float(row,"RAJ2000"),_float(row,"DEJ2000");era,edec,epoch=_float(row,"e_RAJ2000"),_float(row,"e_DEJ2000"),_float(row,"Epoch")
 if None in (ra,dec,era,edec,epoch):raise ValueError("PS1_MEAN_POSITION_OR_EPOCH_INCOMPLETE")
 cov=conservative_unknown_corr_bound_mas2(era,edec,PS1_SYSTEMATIC_FLOOR_MAS)
 return CatalogMeasurement("PANSTARRS_DR1",source_id,ra,dec,epoch,cov.tolist(),"CONSERVATIVE_TRACE_ISOTROPIC_UPPER_BOUND__UNKNOWN_RA_DEC_CORRELATION__20MAS_SYSTEMATIC_FLOOR","EXACT_WEIGHTED_MEAN_EPOCH_MJD_FROM_PS1_EPOCH","PS1_WEIGHTED_MEAN_SINGLE_EPOCH_DETECTION_POSITION_AT_EPOCHMEAN",{"catalog_id":"II/349/ps1","query_mode":row.get("_janus_query_mode"),"raw_epoch_field":"Epoch","raw_epoch_mjd":epoch,"raw_statistical_errors_arcsec":{"ra":era,"dec":edec},"systematic_floor_mas":PS1_SYSTEMATIC_FLOOR_MAS,"unknown_correlation_policy":"TRACE_TIMES_IDENTITY_UPPER_BOUND_NOT_ZERO_CORRELATION"})

def _sdss_sql(sql:str,label:str)->tuple[Dict[str,Any],Dict[str,Any]]:
 resp=requests.get(SDSS_SQL_ENDPOINT,params={"cmd":sql,"format":"csv"},timeout=45);meta={"label":label,"endpoint":SDSS_SQL_ENDPOINT,"request_url":resp.url,"http_status":resp.status_code,"sql":sql};resp.raise_for_status();text=resp.text.strip();lines=[x for x in text.splitlines() if x.strip() and not x.lstrip().startswith("#")];rows=list(csv.DictReader(io.StringIO("\n".join(lines)))) if lines else []
 if not rows:raise RuntimeError(f"SDSS_NO_ROW:{label}:{text[:500]}")
 return rows[0],meta

def query_sdss_exact(objid:str,center:list[float])->Dict[str,Any]:
 oid=int(objid);receipts=[]
 obj,objmeta=_sdss_sql(f"SELECT TOP 1 objID,ra,dec,fieldID,run,rerun,camcol,field FROM PhotoObjAll WHERE objID={oid}","PHOTOOBJ_EXACT_ID");receipts.append(objmeta)
 if str(_get(obj,"objID")).strip()!=str(objid):raise RuntimeError("SDSS_OBJID_MISMATCH")
 fieldid=str(_get(obj,"fieldID")).strip()
 field,fieldmeta=_sdss_sql(f"SELECT TOP 1 fieldID,mjd_u,mjd_g,mjd_r,mjd_i,mjd_z FROM [Field] WHERE fieldID={int(fieldid)}","FIELD_EXACT_FIELDID");receipts.append(fieldmeta)
 vr=query_vizier_exact("SDSS_DR16",objid,center)
 vra,vdec=_float(vr,"RA_ICRS","RAdeg"),_float(vr,"DE_ICRS","DEdeg");era,edec=_float(vr,"e_RA_ICRS","e_RAdeg"),_float(vr,"e_DE_ICRS","e_DEdeg");rmjd=_float(vr,"rMJD");mean_mjd=_float(vr,"MJD");sky_r=_float(field,"mjd_r")
 if None in (vra,vdec,era,edec,rmjd,sky_r):raise ValueError(f"SDSS_VIZIER_POSITION_ERROR_OR_RMJD_INCOMPLETE:{sorted(vr.keys())}")
 if abs(rmjd-sky_r)>0.001:raise ValueError(f"SDSS_RMJD_CROSSCHECK_FAILED:VIZIER={rmjd}:SKYSERVER={sky_r}")
 if str(_get(vr,"fieldID")).strip() not in ("",fieldid):raise ValueError("SDSS_FIELDID_CROSSCHECK_FAILED")
 out={**obj,**field,"RA_catalog":vra,"DE_catalog":vdec,"e_RA_catalog":era,"e_DE_catalog":edec,"rMJD_catalog":rmjd,"MJD_catalog_mean":mean_mjd,"_janus_query_mode":"EXACT_OBJID_SKYSERVER_FIELD_PLUS_EXACT_OBJID_VIZIER_DR16","_janus_query_receipts":receipts,"_janus_vizier_query_mode":vr.get("_janus_query_mode")};return out

def sdss_measurement(source_id:str,row:Dict[str,Any])->CatalogMeasurement:
 ra,dec=_float(row,"RA_catalog"),_float(row,"DE_catalog");era,edec=_float(row,"e_RA_catalog"),_float(row,"e_DE_catalog");mjdr=_float(row,"rMJD_catalog")
 if None in (ra,dec,era,edec,mjdr):raise ValueError("SDSS_RBAND_POSITION_EPOCH_OR_ERRORS_INCOMPLETE")
 cov=conservative_unknown_corr_bound_mas2(era,edec,0.0)
 return CatalogMeasurement("SDSS_DR16",source_id,ra,dec,mjdr,cov.tolist(),"CONSERVATIVE_TRACE_ISOTROPIC_UPPER_BOUND__UNKNOWN_SDSS_RA_DEC_CORRELATION","EXACT_R_BAND_FIELD_EPOCH_MJD_R_CROSSCHECKED_SKYSERVER_VIZIER","SDSS_DR16_CATALOG_POSITION_WITH_R_BAND_FIELD_EPOCH",{"catalog_id":"V/154/sdss16","query_mode":row.get("_janus_query_mode"),"vizier_query_mode":row.get("_janus_vizier_query_mode"),"field_id":_get(row,"fieldID"),"band_mjd":{b:_float(row,f"mjd_{b}") for b in "ugriz"},"vizier_rMJD":mjdr,"vizier_mean_MJD":_float(row,"MJD_catalog_mean"),"raw_statistical_errors_arcsec":{"ra":era,"dec":edec},"unknown_correlation_policy":"TRACE_TIMES_IDENTITY_UPPER_BOUND_NOT_ZERO_CORRELATION","query_receipts":row.get("_janus_query_receipts",[])})

def _resolve_measurement(catalog:str,source_id:str,center:list[float])->CatalogMeasurement:
 if catalog=="SDSS_DR16":return sdss_measurement(source_id,query_sdss_exact(source_id,center))
 row=query_vizier_exact(catalog,source_id,center)
 return allwise_measurement(source_id,row) if catalog=="ALLWISE" else twomass_measurement(source_id,row) if catalog=="2MASS_PSC" else ps1_measurement(source_id,row) if catalog=="PANSTARRS_DR1" else (_ for _ in ()).throw(KeyError(catalog))

def run()->Dict[str,Any]:
 groups={};resolved=0;required=sum(len(g["catalogs"]) for g in GROUPS.values())
 for group_name,cfg in GROUPS.items():
  try:state=state_from_features(cfg["gaia_source_id"],query_gaia_live(cfg["gaia_source_id"]))
  except Exception as exc:groups[group_name]={"status":"I_DO_NOT_KNOW","reason":f"GAIA_QUERY_FAILED:{type(exc).__name__}:{exc}"};continue
  measurements=[];audit={}
  for catalog,source_id in cfg["catalogs"].items():
   try:
    m=_resolve_measurement(catalog,source_id,cfg["center"]);measurements.append(m);resolved+=1;audit[catalog]={"source_id":source_id,"status":"EPOCH_AND_POSITION_COVARIANCE_RESOLVED","epoch_mjd":m.epoch_mjd,"epoch_jyear":m.epoch_jyear,"epoch_status":m.epoch_status,"covariance_status":m.covariance_status,"position_semantics":m.position_semantics}
   except Exception as exc:audit[catalog]={"source_id":source_id,"status":"BLOCKED","reason":f"{type(exc).__name__}:{exc}"}
  groups[group_name]={"status":"EVALUATED" if measurements else "I_DO_NOT_KNOW","target":cfg["target"],"gaia_source_id":cfg["gaia_source_id"],"catalog_epoch_audit":audit,"identity_test":evaluate_identity(state,measurements) if measurements else {"status":"I_DO_NOT_KNOW","reason":"NO_EXTERNAL_MEASUREMENTS_RESOLVED"}}
 return {"schema":SCHEMA,"experiment":"LOVE_EDEM_CATALOG_EPOCH_IDENTITY_GEN3","formula":"RESPICIENS_ET_PROSPICIENS_GEN3","run_time_utc":datetime.now(timezone.utc).isoformat(),"groups":groups,"summary":{"required_external_catalog_measurements":required,"resolved_external_catalog_measurements":resolved,"all_required_epochs_resolved":resolved==required,"release_year_substitutions":0,"simulation_count_increased":False},"archive_search_bridge":{"lineage":"TOPA tools/topa_arxiv_gateway.py","role":"DISCOVERY_AND_PROVENANCE_INDEX_ONLY","catalog_values_must_come_from_primary_catalog_or_archive":True,"search_rank_is_truth":False},"epistemic_firewall":{"catalog_agreement_is_love_edem_identity":False,"catalog_agreement_is_anomaly":False,"missing_epoch_may_use_release_year":False,"missing_covariance_may_be_zeroed":False,"mahalanobis_pass_is_identity_proof":False,"negative_result_is_valid":True},"claim_ceiling":"CROSS_CATALOG_COMPATIBILITY_OF_ALREADY_GROUPED_STELLAR_DETECTIONS_ONLY__NOT_LOVE_EDEM_IDENTITY"}

def self_test()->None:
 a=allwise_measurement("W",{"RA_pm":10.0,"DE_pm":20.0,"e_RA_pm":.1,"e_DE_pm":.2,"cosig_pm":.01});t=twomass_measurement("T",{"RAJ2000":10.0,"DEJ2000":20.0,"errMaj":.2,"errMin":.1,"errPA":30.0,"JD":2451545.0});p=ps1_measurement("P",{"RAJ2000":10.0,"DEJ2000":20.0,"e_RAJ2000":.01,"e_DEJ2000":.02,"Epoch":56000.0});assert a.epoch_mjd==55400.0 and abs(t.epoch_mjd-51544.5)<1e-9 and p.epoch_mjd==56000.0 and "CONSERVATIVE" in p.covariance_status;print("LOVE_EDEM_CATALOG_EPOCH_IDENTITY_SELF_TEST=PASS")

def main()->int:
 p=argparse.ArgumentParser();p.add_argument("--output");p.add_argument("--self-test",action="store_true");a=p.parse_args()
 if a.self_test:self_test();return 0
 out=run();text=json.dumps(out,ensure_ascii=False,indent=2,sort_keys=True);Path(a.output).write_text(text+"\n",encoding="utf-8") if a.output else print(text);return 0
if __name__=="__main__":raise SystemExit(main())
