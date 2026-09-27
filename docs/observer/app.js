(() => {
  const RAW = 'https://raw.githubusercontent.com/Hawkar-usls/Janus-Cosmos/main/';
  const UPSTREAM_PLATES = 'https://raw.githubusercontent.com/jannefi/poss1-plate-slice/main/data/plate_manifest.csv';

  const EARTH_LEDGER = RAW + 'domains/earth/NUC_EVENT_LEDGER_v0.1.json';
  const SPACE_LEDGER = RAW + 'domains/space/PALOMAR_TRANSIENT_LEDGER_v0.1.json';
  const BRIDGE_LEDGER = RAW + 'domains/bridge/EARTH_SPACE_LEDGER_JOIN_CONTRACT_v0.1.json';

  const siteCoordinates = {
    'Pantex Plant': { lat: 35.321, lng: -101.563, precision: 'approximate public facility centroid' },
    'Pilgrim Nuclear Power Station': { lat: 41.944, lng: -70.579, precision: 'approximate public facility centroid' },
    'Los Alamos': { lat: 35.881, lng: -106.299, precision: 'approximate city/laboratory area' },
    'Palomar Observatory': { lat: 33.3563, lng: -116.8650, precision: 'observatory site' }
  };

  const fallbackEarthRows = [
    {
      row_id: 'EARTH-EVENT-001-PANTEX-2015-09-01', row_kind: 'EVENT', date: '2015-09-01', site: 'Pantex Plant',
      region: 'Texas, USA', source_grade: 'P0_PRIMARY_OFFICIAL_INCIDENT_RECORD', classification_state: 'UNRESOLVED_SOURCE_IDENTITY',
      witness_channels: ['ground_surveillance_radar', 'protective_force_visual', 'GSR_imagery'],
      source_url: 'https://www.energy.gov/sites/default/files/2026-09/5%20-%20%28UCNI%29%20NNSA-2015-009581%20-%20COR-NPO-20%20SS-992015-642025.pdf'
    },
    {
      row_id: 'EARTH-EVENT-002-PILGRIM-2015-10-01', row_kind: 'EVENT', date: '2015-10-01', site: 'Pilgrim Nuclear Power Station',
      region: 'Massachusetts, USA', source_grade: 'P0_NARA_SERIES_BOUND__DOCUMENT_TEXT_PENDING_PRIMARY_BYTE_BIND',
      classification_state: 'UNRESOLVED_OBSERVATION__FORMAL_REPORTING_THRESHOLD_NOT_MET', witness_channels: ['human_observation', 'law_enforcement_followup']
    },
    {
      row_id: 'EARTH-HIST-001-LOS-ALAMOS-1949', row_kind: 'HISTORICAL_ARCHIVAL_CASE_FAMILY', date: '1949', site: 'Los Alamos',
      region: 'New Mexico, USA', source_grade: 'P0_OFFICIAL_ARCHIVAL_RELEASE', classification_state: 'INSTITUTIONAL_INVESTIGATION__NO_UNIQUE_ORIGIN_ESTABLISHED'
    }
  ];

  const state = {
    earth: null,
    space: null,
    bridge: null,
    earthPoints: [],
    platePoints: [],
    earthGlobe: null,
    spaceGlobe: null,
    mode: 'earth'
  };

  const $ = (id) => document.getElementById(id);

  function tickClock() {
    const now = new Date();
    $('clock').textContent = now.toISOString().replace('T', ' ').slice(0, 19) + ' UTC';
  }
  tickClock();
  setInterval(tickClock, 1000);

  function safeText(v, fallback = '—') {
    return v === undefined || v === null || v === '' ? fallback : String(v);
  }

  function gradeShort(grade) {
    if (!grade) return 'OPEN';
    if (grade.startsWith('P0_')) return 'P0';
    if (grade.includes('PRIMARY')) return 'PRIMARY';
    return grade.slice(0, 12);
  }

  function classificationColor(row) {
    const s = (row.classification_state || '').toUpperCase();
    if (s.includes('POLICY')) return '#ffbd5a';
    if (s.includes('UAS')) return '#9a8cff';
    if (s.includes('UNRESOLVED')) return '#53d7ff';
    return '#65b7ff';
  }

  function rowToPoint(row) {
    const c = siteCoordinates[row.site];
    if (!c) return null;
    return {
      ...row,
      lat: c.lat,
      lng: c.lng,
      coordinate_precision: c.precision,
      color: classificationColor(row),
      size: row.row_kind === 'EVENT' ? 0.34 : 0.26
    };
  }

  async function fetchJSON(url) {
    const r = await fetch(url, { cache: 'no-store' });
    if (!r.ok) throw new Error(`${r.status} ${url}`);
    return r.json();
  }

  async function fetchText(url) {
    const r = await fetch(url, { cache: 'no-store' });
    if (!r.ok) throw new Error(`${r.status} ${url}`);
    return r.text();
  }

  function parsePlateManifest(csv) {
    const lines = csv.trim().split(/\r?\n/);
    const header = lines.shift().split(',');
    const idx = Object.fromEntries(header.map((h, i) => [h.trim(), i]));
    return lines.map(line => {
      const c = line.split(',');
      const ra = Number(c[idx.ra_deg]);
      const dec = Number(c[idx.dec_deg]);
      return {
        plate_id: c[idx.plate_id],
        ra,
        dec,
        lat: dec,
        lng: ra > 180 ? ra - 360 : ra,
        color: dec > 60 ? '#64d8ff' : dec > 20 ? '#3da7ff' : '#8d77ff'
      };
    }).filter(p => Number.isFinite(p.ra) && Number.isFinite(p.dec));
  }

  function pointLabel(p) {
    return `<div class="globe-label"><b>${safeText(p.site)}</b><small>${safeText(p.date || p.row_kind)}</small></div>`;
  }

  function plateLabel(p) {
    return `<div class="celestial-label"><b>${p.plate_id}</b><br>RA ${p.ra.toFixed(2)}° · Dec ${p.dec.toFixed(2)}°</div>`;
  }

  function showDetail(p, domain = 'earth') {
    if (!p) return;
    $('detail-grade').textContent = domain === 'earth' ? gradeShort(p.source_grade) : 'POSS-I';
    if (domain === 'space') {
      $('detail').className = 'detail-content';
      $('detail').innerHTML = `
        <span class="sub">SPACE · SURVEY POSITION</span>
        <h3>${p.plate_id}</h3>
        <div class="detail-table">
          <div><span>RIGHT ASC.</span><b>${p.ra.toFixed(6)}°</b></div>
          <div><span>DECLINATION</span><b>${p.dec.toFixed(6)}°</b></div>
          <div><span>REPRESENTATION</span><b>POSS-I plate center</b></div>
          <div><span>CLAIM CEILING</span><b>Survey position only — not an anomaly claim</b></div>
        </div>
        <a class="source-link" target="_blank" rel="noreferrer" href="https://github.com/jannefi/poss1-plate-slice/blob/main/data/plate_manifest.csv">OPEN PLATE MANIFEST ↗</a>`;
      return;
    }

    const channels = Array.isArray(p.witness_channels) ? p.witness_channels.join(' · ') : safeText(p.row_kind);
    $('detail').className = 'detail-content';
    $('detail').innerHTML = `
      <span class="sub">EARTH · ${safeText(p.row_kind)}</span>
      <h3>${safeText(p.site)}</h3>
      <div class="detail-table">
        <div><span>DATE</span><b>${safeText(p.date || p.interval || p.coverage_start)}</b></div>
        <div><span>REGION</span><b>${safeText(p.region)}</b></div>
        <div><span>STATE</span><b>${safeText(p.classification_state)}</b></div>
        <div><span>WITNESSES</span><b>${channels}</b></div>
        <div><span>COORDS</span><b>${p.lat.toFixed(3)}, ${p.lng.toFixed(3)} · ${p.coordinate_precision}</b></div>
        <div><span>ROW ID</span><b>${safeText(p.row_id)}</b></div>
      </div>
      ${p.source_url ? `<a class="source-link" target="_blank" rel="noreferrer" href="${p.source_url}">OPEN PRIMARY / SOURCE PAGE ↗</a>` : ''}`;
  }

  function buildEventList(points) {
    const list = $('earth-list');
    list.classList.remove('skeleton-list');
    list.innerHTML = '';
    points.forEach(p => {
      const el = document.createElement('div');
      el.className = 'event-card';
      const provisional = (p.source_grade || '').includes('PENDING') || (p.source_grade || '').includes('SERIES_BOUND');
      el.innerHTML = `<div class="row"><b>${safeText(p.site)}</b><small class="badge ${provisional ? 'warn' : ''}">${gradeShort(p.source_grade)}</small></div>
        <div class="row kind"><span>${safeText(p.date || p.row_kind)}</span><small>${safeText(p.region)}</small></div>`;
      el.addEventListener('click', () => {
        showDetail(p, 'earth');
        if (state.earthGlobe) state.earthGlobe.pointOfView({ lat: p.lat, lng: p.lng, altitude: 1.55 }, 900);
      });
      list.appendChild(el);
    });
  }

  function initEarthGlobe(points) {
    const node = $('earth-globe');
    const globe = Globe()(node)
      .backgroundColor('rgba(0,0,0,0)')
      .globeImageUrl('https://unpkg.com/three-globe/example/img/earth-night.jpg')
      .bumpImageUrl('https://unpkg.com/three-globe/example/img/earth-topology.png')
      .showAtmosphere(true)
      .atmosphereColor('#3caeff')
      .atmosphereAltitude(0.18)
      .pointsData(points)
      .pointLat('lat').pointLng('lng')
      .pointColor('color')
      .pointAltitude(0.012)
      .pointRadius('size')
      .pointResolution(16)
      .htmlElementsData(points)
      .htmlLat('lat').htmlLng('lng').htmlAltitude(0.055)
      .htmlElement(p => {
        const el = document.createElement('div');
        el.innerHTML = pointLabel(p);
        el.style.pointerEvents = 'auto';
        el.style.cursor = 'pointer';
        el.onclick = () => showDetail(p, 'earth');
        return el;
      })
      .onPointClick(p => showDetail(p, 'earth'))
      .onPointHover(p => { if (p) showDetail(p, 'earth'); })
      .pointOfView({ lat: 33, lng: -55, altitude: 1.9 });

    globe.controls().autoRotate = true;
    globe.controls().autoRotateSpeed = 0.28;
    globe.controls().enableDamping = true;
    state.earthGlobe = globe;
  }

  function initSpaceGlobe(plates) {
    const node = $('space-globe');
    const globe = Globe()(node)
      .backgroundColor('rgba(0,0,0,0)')
      .showGlobe(true)
      .showAtmosphere(true)
      .atmosphereColor('#604dff')
      .atmosphereAltitude(0.12)
      .pointsData(plates)
      .pointLat('lat').pointLng('lng')
      .pointColor('color')
      .pointAltitude(0.026)
      .pointRadius(0.12)
      .pointResolution(8)
      .pointLabel(plateLabel)
      .onPointClick(p => showDetail(p, 'space'))
      .onPointHover(p => { if (p) showDetail(p, 'space'); })
      .pointOfView({ lat: 30, lng: 25, altitude: 1.72 });

    const mat = globe.globeMaterial();
    if (mat) {
      mat.color = new THREE.Color('#091326');
      mat.opacity = 0.28;
      mat.transparent = true;
      mat.emissive = new THREE.Color('#050b1c');
      mat.emissiveIntensity = 0.6;
    }
    globe.controls().autoRotate = true;
    globe.controls().autoRotateSpeed = -0.2;
    globe.controls().enableDamping = true;
    state.spaceGlobe = globe;
  }

  function resizeGlobes() {
    const stage = document.querySelector('.stage');
    if (!stage) return;
    const w = stage.clientWidth;
    const h = stage.clientHeight;
    if (state.earthGlobe) state.earthGlobe.width(w).height(h);
    if (state.spaceGlobe) state.spaceGlobe.width(w).height(h);
  }

  function setMode(mode) {
    state.mode = mode;
    document.querySelectorAll('.mode').forEach(b => b.classList.toggle('active', b.dataset.mode === mode));
    $('earth-globe').classList.toggle('active', mode === 'earth');
    $('space-globe').classList.toggle('active', mode === 'space');
    $('bridge-view').classList.toggle('active', mode === 'bridge');

    if (mode === 'earth') {
      $('stage-eyebrow').textContent = 'TERRESTRIAL WITNESS LANE';
      $('stage-title').textContent = 'Interactive Earth';
      $('stage-big').textContent = state.earthPoints.length || '—';
      $('stage-small').textContent = 'mapped source-bound observation sites';
      $('hud-source').textContent = 'EARTH LEDGER';
    } else if (mode === 'space') {
      $('stage-eyebrow').textContent = 'CELESTIAL SURVEY LANE';
      $('stage-title').textContent = 'POSS-I Observation Sphere';
      $('stage-big').textContent = state.platePoints.length || '—';
      $('stage-small').textContent = 'public plate centers';
      $('hud-source').textContent = 'SPACE LEDGER';
    } else {
      $('stage-eyebrow').textContent = 'CONTROLLED CROSS-DOMAIN LANE';
      $('stage-title').textContent = 'Earth ↔ Space Bridge';
      $('stage-big').textContent = 'L2';
      $('stage-small').textContent = 'exact-day join gate';
      $('hud-source').textContent = 'JOIN CONTRACT';
    }
    setTimeout(resizeGlobes, 80);
  }

  document.querySelectorAll('.mode').forEach(btn => btn.addEventListener('click', () => setMode(btn.dataset.mode)));
  window.addEventListener('resize', () => requestAnimationFrame(resizeGlobes));

  async function load() {
    $('load-state').textContent = 'LOADING';
    const tasks = await Promise.allSettled([
      fetchJSON(EARTH_LEDGER),
      fetchJSON(SPACE_LEDGER),
      fetchJSON(BRIDGE_LEDGER),
      fetchText(UPSTREAM_PLATES)
    ]);

    state.earth = tasks[0].status === 'fulfilled' ? tasks[0].value : { rows: fallbackEarthRows };
    state.space = tasks[1].status === 'fulfilled' ? tasks[1].value : null;
    state.bridge = tasks[2].status === 'fulfilled' ? tasks[2].value : null;

    const earthRows = Array.isArray(state.earth?.rows) ? state.earth.rows : fallbackEarthRows;
    state.earthPoints = earthRows.map(rowToPoint).filter(Boolean);

    if (tasks[3].status === 'fulfilled') {
      state.platePoints = parsePlateManifest(tasks[3].value);
    }

    $('earth-count').textContent = `${state.earthPoints.length} MAPPED`;
    buildEventList(state.earthPoints);

    const cohort = state.space?.rows?.find(r => r.row_kind === 'ANALYSIS_COHORT');
    const dataset = state.space?.rows?.find(r => r.row_kind === 'PUBLIC_DATASET_RELEASE');
    $('m-plates').textContent = safeText(dataset?.plate_count || cohort?.plates_in_window);
    $('m-tiles').textContent = safeText(dataset?.tile_rows);
    $('m-candidates').textContent = Number(dataset?.catalogue_rows || cohort?.candidate_rows_on_in_window_plates || 0).toLocaleString() || '—';
    $('m-days').textContent = safeText(cohort?.unique_observation_dates);
    $('space-count').textContent = state.platePoints.length ? `${state.platePoints.length} SKY POSITIONS` : 'SOURCE OPEN';

    if (state.bridge) {
      const l2 = state.bridge.join_levels?.find(x => x.level === 'L2_EXACT_DAY');
      $('bridge-status').textContent = l2?.allowed ? 'L2 EXACT-DAY JOIN: OPEN' : 'L2 EXACT-DAY JOIN: BLOCKED';
      $('bridge-reason').textContent = l2?.current_blocker || l2?.unlock || 'Bridge contract loaded.';
    }

    initEarthGlobe(state.earthPoints);
    if (state.platePoints.length) initSpaceGlobe(state.platePoints);
    resizeGlobes();
    setMode('earth');

    $('load-state').textContent = tasks.every(t => t.status === 'fulfilled') ? 'LIVE' : 'PARTIAL';
    $('load-state').style.color = tasks.every(t => t.status === 'fulfilled') ? '#48e59a' : '#ffbd5a';
  }

  load().catch(err => {
    console.error(err);
    $('load-state').textContent = 'DEGRADED';
    $('load-state').style.color = '#ff6b6b';
  });
})();
