#!/usr/bin/env python3
"""Read-only evidence acquisition for one YouTube video, no content claims without captions."""
from pathlib import Path
import hashlib, json, os, re, sys, traceback
from urllib.parse import quote
import requests

VIDEO_ID = "oMPHmMnS0Mc"
URL = "https://www.youtube.com/watch?v=" + VIDEO_ID
OUT = Path("data/cousteau/video_intake") / f"youtube-{VIDEO_ID}-run-{os.environ.get('GITHUB_RUN_ID','local')}.json"
HEAD = {"User-Agent":"Mozilla/5.0 (JANUS Cousteau video-source verification)"}

def trim(txt,limit):
    return re.sub(r"[\\r\\x00]", "",str(txt))[:limit]

def safe_get(url):
    try:
        r=requests.get(url,headers=HEAD,timeout=12)
        return {"http_status":r.status_code,"text":r.text[:120000]}
    except Exception as e:
        return {"error":type(e).__name__+":"+trim(e,120)}

receipt={"schema":"janus.cousteau.youtube_evidence_intake.v1",
 "video_id":VIDEO_ID,"source_url":URL,"state":"EVIDENCE_PENDING",
 "claims_permitted":False,"content_analysis_complete":False,"metadata_sources":[],
 "captions_found":False,"visual_frame_inspection":False,"report":[]}

try:
    oembed_url="https://www.youtube.com/oembed?url="+quote(URL,safe="")+"&format=json"
    for provider,u in [("youtube_oembed",oembed_url),
       ("noembed","https://noembed.com/embed?url="+quote(URL,safe=""))]:
        r=safe_get(u)
        row={"provider":provider,"http_status":r.get("http_status"),"error":r.get("error")}
        if r.get("http_status")==200:
            try:
                obj=json.loads(r["text"])
                row.update({"title":trim(obj.get("title"),300),"author_name":trim(obj.get("author_name"),200)})
            except Exception as e:row["error"]="json_parse:"+type(e).__name__
        receipt["metadata_sources"].append(row)
    try:
        import yt_dlp
        opts={"skip_download":True,"quiet":True,"no_warnings":True,"ignoreerrors":False,
          "socket_timeout":12,"retries":1,"extractor_retries":1,"noplaylist":True}
        with yt_dlp.YoutubeDL(opts) as ydl:
            inf=ydl.extract_info(URL,download=False)
            if inf:
                receipt["yt_dlp"]={
                  "title":trim(inf.get("title"),300),"channel":trim(inf.get("channel"),200),
                  "description":trim(inf.get("description"),3500),
                  "duration_seconds":inf.get("duration"),"upload_date":inf.get("upload_date"),
                  "chapters":[{"start_time":c.get("start_time"),"title":trim(c.get("title"),200)} for c in (inf.get("chapters") or [])[:30]],
                  "subtitles_languages":list((inf.get("subtitles") or {}).keys()),
                  "automatic_captions_languages":list((inf.get("automatic_captions") or {}).keys())[:30]}
    except Exception as e:
        receipt["yt_dlp_error"]=type(e).__name__+":"+trim(e,450)
    segments=[]
    try:
        from youtube_transcript_api import YouTubeTranscriptApi
        cli=YouTubeTranscriptApi()
        listing=cli.list(VIDEO_ID)
        available=list(listing)
        order=["ru","en"]
        for tr in sorted(available,key=lambda z:(0 if z.language_code in order else 1,z.is_generated)):
            try:
                track=tr.fetch()
                for seg in track:
                    segments.append({"start":round(float(seg.start),2),"duration":round(float(seg.duration),2),"text":trim(seg.text,500)})
                receipt["transcript_language"]=tr.language_code
                receipt["transcript_generated"]=tr.is_generated
                break
            except Exception as e:
                receipt.setdefault("caption_errors",[]).append(type(e).__name__+":"+trim(e,170))
    except Exception as e:
        receipt["caption_list_error"]=type(e).__name__+":"+trim(e,400)
    receipt["captions_found"]=bool(segments)
    if segments:
        receipt["caption_count"]=len(segments)
        receipt["caption_sha256"]=hashlib.sha256(json.dumps(segments,sort_keys=True).encode()).hexdigest()
        # Representative time-distributed extracts for evidence inspection, not the entire transcript.
        indices=sorted(set(int(i*(len(segments)-1)/min(len(segments)-1,31)) for i in range(min(len(segments),32))))
        receipt["caption_sample"]=[segments[i] for i in indices]
        receipt["content_analysis_complete"]=False
        receipt["state"]="CAPTIONS_RETRIEVED__PENDING_HUMAN_OR_MODEL_REVIEW"
    elif receipt.get("yt_dlp",{}).get("title") or any(x.get("title") for x in receipt["metadata_sources"]):
        receipt["state"]="METADATA_ONLY__NO_VERIFIED_CAPTIONS"
    else:receipt["state"]="SOURCE_UNAVAILABLE_NO_CLAIMS"
except Exception as e:
    receipt["state"]="ERROR_NO_CLAIMS"
    receipt["error"]=type(e).__name__+":"+trim(e,700)
OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(receipt,ensure_ascii=False,indent=2,sort_keys=True)+"\\n",encoding="utf-8")
print("COUSTEAU_VIDEO_EVIDENCE_BEGIN")
print(json.dumps(receipt,ensure_ascii=False,sort_keys=True))
print("COUSTEAU_VIDEO_EVIDENCE_END")
print("EVIDENCE_SAVED",str(OUT),receipt["state"])
