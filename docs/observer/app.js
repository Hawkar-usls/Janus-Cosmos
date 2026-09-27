(() => {
  'use strict';

  const DATA = {
    earth: './data/earth-ledger.json',
    space: './data/space-ledger.json',
    bridge: './data/bridge-ledger.json',
    plates: './data/plate_manifest.csv'
  };

  const siteCoordinates = {
    'Pantex Plant': { lat: 35.321, lng: -101.563, precision: 'approximate public facility centroid' },
    'Pilgrim Nuclear Power Station': { lat: 41.944, lng: -70.579, precision: 'approximate public facility centroid' },
    'Los Alamos': { lat: 35.881, lng: -106.299, precision: 'approximate city/laboratory area' },
    'Palomar Observatory': { lat: 33.3563, lng: -116.8650, precision: 'observatory site' }
  };

  const fallbackEarthRows = [
    {
      row_id: 'EARTH-EVENT-001-PANTEX-2015-09-01', row_kind: 'EVENT', date: '2015-09-01', site: 'Pantex Plant', region: 'Texas, USA',
      source_grade: 'P0_PRIMARY_OFFICIAL_INCIDENT_RECORD', classification_state: 'UNRESOLVED_SOURCE_IDENTITY',
      witness_channels: ['ground_surveillance_radar', 'protective_force_visual', 'GSR_imagery'],
      source_url: 'https://www.energy.gov/sites/default/files/2026-09/5%20-%20%28UCNI%29%20NNSA-2015-009581%20-%20COR-NPO-20%20SS-992015-642025.pdf'
    },
    {
      row_id: 'EARTH-EVENT-002-PILGRIM-2015-10-01', row_kind: 'EVENT', date: '2015-10-01', site: 'Pilgrim Nuclear Power Station', region: 'Massachusetts, USA',
      source_grade: 'P0_NARA_SERIES_BOUND__DOCUMENT_TEXT_PENDING_PRIMARY_BYTE_BIND',
      classification_state: 'UNRESOLVED_OBSERVATION__FORMAL_REPORTING_THRESHOLD_NOT_MET',
      witness_channels: ['human_observation', 'law_enforcement_followup']
    },
    {
      row_id: 'EARTH-HIST-001-LOS-ALAMOS-1949', row_kind: 'HISTORICAL_ARCHIVAL_CASE_FAMILY', date: '1949', site: 'Los Alamos', region: 'New Mexico, USA',
      source_grade: 'P0_OFFICIAL_ARCHIVAL_RELEASE', classification_state: 'INSTITUTIONAL_INVESTIGATION__NO_UNIQUE_ORIGIN_ESTABLISHED'
    }
  ];

  const fallbackPlateCSV = `plate_id,ra_deg,dec_deg
XE002,2.348777,84.754856
XE003,44.912000,84.549934
XE004,83.899465,84.044706
XE005,120.975302,83.586154
XE006,154.842758,83.263035
XE007,188.661150,83.255440
XE008,221.640565,83.384269`;

  const continents = [
    [[72,-165],[62,-150],[53,-135],[48,-125],[35,-120],[23,-108],[19,-98],[25,-82],[42,-70],[52,-58],[61,-74],[70,-105],[72,-165]],
    [[12,-80],[4,-76],[-12,-72],[-30,-70],[-55,-68],[-43,-55],[-22,-45],[0,-50],[10,-65],[12,-80]],
    [[70,-10],[61,17],[49,30],[38,24],[35,8],[35,-10],[19,-17],[2,-8],[-18,12],[-35,20],[-34,35],[-12,45],[11,50],[31,35],[43,25],[52,10],[70,-10]],
    [[70,28],[72,80],[62,128],[52,155],[36,142],[22,112],[8,80],[24,57],[39,45],[55,48],[70,28]],
    [[-10,112],[-21,114],[-35,116],[-42,145],[-25,154],[-11,145],[-10,112]],
    [[83,-45],[75,-18],[66,-35],[62,-52],[70,-65],[83,-45]]
  ];

  const state = {
    earth: null,
    space: null,
    bridge: null,
    earthPoints: [],
    platePoints: [],
    earthSphere: null,
    spaceSphere: null,
    mode: 'earth'
  };

  const $ = id => document.getElementById(id);
  const rad = d => d * Math.PI / 180;
  const clamp = (v, a, b) => Math.max(a, Math.min(b, v));

  function tickClock() {
    $('clock').textContent = new Date().toISOString().replace('T', ' ').slice(0, 19) + ' UTC';
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
      size: row.row_kind === 'EVENT' ? 7 : 5
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
    if (lines.length < 2) return [];
    const header = lines.shift().split(',');
    const idx = Object.fromEntries(header.map((h, i) => [h.trim(), i]));
    return lines.map(line => {
      const c = line.split(',');
      const ra = Number(c[idx.ra_deg]);
      const dec = Number(c[idx.dec_deg]);
      return {
        plate_id: c[idx.plate_id],
        ra, dec,
        lat: dec,
        lng: ra > 180 ? ra - 360 : ra,
        color: dec > 60 ? '#6de3ff' : dec > 20 ? '#3fa8ff' : '#9c7bff',
        size: 2.7
      };
    }).filter(p => p.plate_id && Number.isFinite(p.ra) && Number.isFinite(p.dec));
  }

  function seededStars(count) {
    let seed = 0x1a2b3c4d;
    const rnd = () => {
      seed = (1664525 * seed + 1013904223) >>> 0;
      return seed / 4294967296;
    };
    return Array.from({ length: count }, () => ({
      lat: Math.asin(rnd() * 2 - 1) * 180 / Math.PI,
      lng: rnd() * 360 - 180,
      alpha: .15 + rnd() * .65,
      size: .45 + rnd() * 1.4
    }));
  }

  class JanusSphere {
    constructor(node, options = {}) {
      this.node = node;
      this.mode = options.mode || 'earth';
      this.points = [];
      this.screenPoints = [];
      this.yaw = this.mode === 'earth' ? rad(-40) : rad(15);
      this.pitch = this.mode === 'earth' ? rad(-8) : rad(8);
      this.zoom = 1;
      this.dragging = false;
      this.lastX = 0;
      this.lastY = 0;
      this.lastInteraction = 0;
      this.hovered = null;
      this.selected = null;
      this.onHover = options.onHover || (() => {});
      this.onSelect = options.onSelect || (() => {});
      this.stars = seededStars(340);

      this.canvas = document.createElement('canvas');
      this.canvas.setAttribute('aria-hidden', 'true');
      Object.assign(this.canvas.style, { position: 'absolute', inset: '0', width: '100%', height: '100%', display: 'block', touchAction: 'none' });
      node.replaceChildren(this.canvas);
      this.ctx = this.canvas.getContext('2d', { alpha: true });

      this.bind();
      this.resize();
      this.ro = new ResizeObserver(() => this.resize());
      this.ro.observe(node);
      this.loop = this.loop.bind(this);
      requestAnimationFrame(this.loop);
    }

    setPoints(points) {
      this.points = Array.isArray(points) ? points : [];
      return this;
    }

    resize() {
      const rect = this.node.getBoundingClientRect();
      const dpr = Math.min(window.devicePixelRatio || 1, 2);
      this.w = Math.max(1, rect.width);
      this.h = Math.max(1, rect.height);
      this.canvas.width = Math.floor(this.w * dpr);
      this.canvas.height = Math.floor(this.h * dpr);
      this.canvas.style.width = `${this.w}px`;
      this.canvas.style.height = `${this.h}px`;
      this.ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      this.radius = Math.min(this.w, this.h) * (this.mode === 'earth' ? .37 : .39) * this.zoom;
      this.cx = this.w * .5;
      this.cy = this.h * .54;
    }

    bind() {
      this.canvas.addEventListener('pointerdown', e => {
        this.dragging = true;
        this.lastX = e.clientX;
        this.lastY = e.clientY;
        this.lastInteraction = performance.now();
        this.canvas.setPointerCapture?.(e.pointerId);
      });
      this.canvas.addEventListener('pointerup', e => {
        const moved = Math.hypot(e.clientX - this.lastX, e.clientY - this.lastY);
        this.dragging = false;
        this.lastInteraction = performance.now();
        if (moved < 8) this.pick(e.offsetX, e.offsetY, true);
      });
      this.canvas.addEventListener('pointercancel', () => { this.dragging = false; });
      this.canvas.addEventListener('pointermove', e => {
        if (this.dragging) {
          const dx = e.clientX - this.lastX;
          const dy = e.clientY - this.lastY;
          this.lastX = e.clientX;
          this.lastY = e.clientY;
          this.yaw += dx * .0075;
          this.pitch = clamp(this.pitch + dy * .0065, -1.25, 1.25);
          this.lastInteraction = performance.now();
        } else {
          this.pick(e.offsetX, e.offsetY, false);
        }
      });
      this.canvas.addEventListener('wheel', e => {
        e.preventDefault();
        this.zoom = clamp(this.zoom * (e.deltaY > 0 ? .92 : 1.08), .72, 1.45);
        this.lastInteraction = performance.now();
        this.resize();
      }, { passive: false });
    }

    pick(x, y, select) {
      let best = null;
      let bestD = 18;
      for (const s of this.screenPoints) {
        const d = Math.hypot(x - s.x, y - s.y);
        if (d < bestD) { best = s.point; bestD = d; }
      }
      if (select && best) {
        this.selected = best;
        this.onSelect(best);
      } else if (!select && best !== this.hovered) {
        this.hovered = best;
        if (best) this.onHover(best);
      }
      this.canvas.style.cursor = best ? 'pointer' : (this.dragging ? 'grabbing' : 'grab');
    }

    xyz(lat, lng) {
      const la = rad(lat), lo = rad(lng);
      let x = Math.cos(la) * Math.sin(lo);
      let y = Math.sin(la);
      let z = Math.cos(la) * Math.cos(lo);
      const cy = Math.cos(this.yaw), sy = Math.sin(this.yaw);
      const x1 = x * cy + z * sy;
      const z1 = -x * sy + z * cy;
      x = x1; z = z1;
      const cp = Math.cos(this.pitch), sp = Math.sin(this.pitch);
      const y1 = y * cp - z * sp;
      const z2 = y * sp + z * cp;
      return { x, y: y1, z: z2 };
    }

    project(lat, lng, altitude = 0) {
      const p = this.xyz(lat, lng);
      const s = 1 + altitude;
      return { x: this.cx + p.x * this.radius * s, y: this.cy - p.y * this.radius * s, z: p.z, visible: p.z > -0.02 };
    }

    pathLat(lat, color, alpha = .24) {
      const c = this.ctx;
      let open = false;
      c.beginPath();
      for (let lng = -180; lng <= 180; lng += 3) {
        const p = this.project(lat, lng);
        if (p.z > 0) {
          if (!open) { c.moveTo(p.x, p.y); open = true; } else c.lineTo(p.x, p.y);
        } else open = false;
      }
      c.strokeStyle = color;
      c.globalAlpha = alpha;
      c.lineWidth = 1;
      c.stroke();
      c.globalAlpha = 1;
    }

    pathLng(lng, color, alpha = .22) {
      const c = this.ctx;
      let open = false;
      c.beginPath();
      for (let lat = -88; lat <= 88; lat += 2) {
        const p = this.project(lat, lng);
        if (p.z > 0) {
          if (!open) { c.moveTo(p.x, p.y); open = true; } else c.lineTo(p.x, p.y);
        } else open = false;
      }
      c.strokeStyle = color;
      c.globalAlpha = alpha;
      c.lineWidth = 1;
      c.stroke();
      c.globalAlpha = 1;
    }

    drawContinents() {
      const c = this.ctx;
      c.save();
      c.strokeStyle = '#60dfff';
      c.lineWidth = 1.35;
      c.shadowColor = '#1ca9ff';
      c.shadowBlur = 7;
      c.globalAlpha = .6;
      for (const poly of continents) {
        c.beginPath();
        let open = false;
        for (const [lat, lng] of poly) {
          const p = this.project(lat, lng, .004);
          if (p.z > .02) {
            if (!open) { c.moveTo(p.x, p.y); open = true; } else c.lineTo(p.x, p.y);
          } else open = false;
        }
        c.stroke();
      }
      c.restore();
    }

    drawStars() {
      const c = this.ctx;
      c.save();
      for (const s of this.stars) {
        const p = this.project(s.lat, s.lng, -.015);
        if (p.z <= 0) continue;
        c.globalAlpha = s.alpha * (.45 + p.z * .55);
        c.fillStyle = p.z > .72 ? '#dff7ff' : '#6ea9ff';
        c.beginPath();
        c.arc(p.x, p.y, s.size, 0, Math.PI * 2);
        c.fill();
      }
      c.restore();
    }

    draw() {
      if (!this.ctx || !this.w || !this.h) return;
      const c = this.ctx;
      c.clearRect(0, 0, this.w, this.h);
      const R = this.radius;
      const glow = c.createRadialGradient(this.cx - R * .28, this.cy - R * .28, R * .06, this.cx, this.cy, R * 1.18);
      if (this.mode === 'earth') {
        glow.addColorStop(0, '#1765a5'); glow.addColorStop(.38, '#092b54'); glow.addColorStop(.8, '#041126'); glow.addColorStop(1, '#010611');
      } else {
        glow.addColorStop(0, '#2a245d'); glow.addColorStop(.42, '#111a48'); glow.addColorStop(.82, '#060b23'); glow.addColorStop(1, '#01040e');
      }
      c.save();
      c.shadowColor = this.mode === 'earth' ? '#159cff' : '#745cff';
      c.shadowBlur = 42;
      c.fillStyle = glow;
      c.beginPath(); c.arc(this.cx, this.cy, R, 0, Math.PI * 2); c.fill(); c.restore();

      c.save();
      c.beginPath(); c.arc(this.cx, this.cy, R - 1, 0, Math.PI * 2); c.clip();
      if (this.mode === 'space') this.drawStars();
      const grid = this.mode === 'earth' ? '#45bfff' : '#836fff';
      for (let lat = -60; lat <= 60; lat += 30) this.pathLat(lat, grid, lat === 0 ? .34 : .17);
      for (let lng = -150; lng <= 180; lng += 30) this.pathLng(lng, grid, .15);
      if (this.mode === 'earth') this.drawContinents();
      const limb = c.createLinearGradient(this.cx - R, this.cy, this.cx + R, this.cy);
      limb.addColorStop(0, 'rgba(0,0,0,.55)'); limb.addColorStop(.55, 'rgba(0,0,0,0)'); limb.addColorStop(1, 'rgba(255,255,255,.05)');
      c.fillStyle = limb; c.fillRect(this.cx - R, this.cy - R, R * 2, R * 2); c.restore();

      c.save();
      c.strokeStyle = this.mode === 'earth' ? 'rgba(83,215,255,.72)' : 'rgba(139,113,255,.74)';
      c.lineWidth = 1.5; c.shadowColor = this.mode === 'earth' ? '#42d0ff' : '#715cff'; c.shadowBlur = 18;
      c.beginPath(); c.arc(this.cx, this.cy, R, 0, Math.PI * 2); c.stroke(); c.restore();

      this.screenPoints = [];
      const visible = this.points.map(point => ({ point, p: this.project(point.lat, point.lng, .025) }))
        .filter(o => o.p.z > 0).sort((a, b) => a.p.z - b.p.z);
      for (const { point, p } of visible) {
        const depth = .55 + p.z * .75;
        const base = point.size || (this.mode === 'earth' ? 5 : 2.5);
        const size = base * depth;
        c.save(); c.globalAlpha = .65 + p.z * .35; c.fillStyle = point.color || '#53d7ff'; c.shadowColor = point.color || '#53d7ff'; c.shadowBlur = this.mode === 'earth' ? 18 : 10;
        c.beginPath(); c.arc(p.x, p.y, size, 0, Math.PI * 2); c.fill();
        c.globalAlpha = .78; c.strokeStyle = '#ffffff'; c.lineWidth = .8; c.beginPath(); c.arc(p.x, p.y, size + 3, 0, Math.PI * 2); c.stroke(); c.restore();
        this.screenPoints.push({ point, x: p.x, y: p.y });
        if (this.mode === 'earth' && (point.row_kind === 'EVENT' || point === this.selected)) {
          c.save(); c.fillStyle = '#e8f8ff'; c.font = '600 11px system-ui, sans-serif'; c.shadowColor = '#00101f'; c.shadowBlur = 5;
          c.fillText(point.site, p.x + 11, p.y - 7); c.fillStyle = '#79bde4'; c.font = '9px system-ui, sans-serif'; c.fillText(safeText(point.date || point.row_kind), p.x + 11, p.y + 6); c.restore();
        }
      }
      if (this.hovered && this.mode === 'space') {
        const s = this.screenPoints.find(x => x.point === this.hovered);
        if (s) {
          c.save(); c.font = '600 10px system-ui, sans-serif';
          const label = `${this.hovered.plate_id} · RA ${this.hovered.ra.toFixed(2)}° · Dec ${this.hovered.dec.toFixed(2)}°`;
          const tw = c.measureText(label).width;
          c.fillStyle = 'rgba(3,10,24,.88)'; c.fillRect(s.x + 10, s.y - 26, tw + 16, 22);
          c.strokeStyle = 'rgba(83,215,255,.35)'; c.strokeRect(s.x + 10, s.y - 26, tw + 16, 22);
          c.fillStyle = '#dff6ff'; c.fillText(label, s.x + 18, s.y - 11); c.restore();
        }
      }
    }

    loop(t) {
      if (!this.dragging && t - this.lastInteraction > 1800) this.yaw += this.mode === 'earth' ? .0006 : -.00042;
      this.draw();
      requestAnimationFrame(this.loop);
    }

    focus(point) {
      if (!point) return;
      this.selected = point;
      this.yaw = -rad(point.lng);
      this.pitch = rad(point.lat) * .42;
      this.lastInteraction = performance.now();
    }
  }

  function showDetail(p, domain = 'earth') {
    if (!p) return;
    $('detail-grade').textContent = domain === 'earth' ? gradeShort(p.source_grade) : 'POSS-I';
    if (domain === 'space') {
      $('detail').className = 'detail-content';
      $('detail').innerHTML = `<span class="sub">SPACE · SURVEY POSITION</span><h3>${safeText(p.plate_id)}</h3><div class="detail-table"><div><span>RIGHT ASC.</span><b>${Number(p.ra).toFixed(6)}°</b></div><div><span>DECLINATION</span><b>${Number(p.dec).toFixed(6)}°</b></div><div><span>REPRESENTATION</span><b>POSS-I plate center</b></div><div><span>CLAIM CEILING</span><b>Survey position only — not an anomaly claim</b></div></div><a class="source-link" target="_blank" rel="noreferrer" href="https://github.com/jannefi/poss1-plate-slice/blob/main/data/plate_manifest.csv">OPEN PLATE MANIFEST ↗</a>`;
      return;
    }
    const channels = Array.isArray(p.witness_channels) ? p.witness_channels.join(' · ') : safeText(p.row_kind);
    $('detail').className = 'detail-content';
    $('detail').innerHTML = `<span class="sub">EARTH · ${safeText(p.row_kind)}</span><h3>${safeText(p.site)}</h3><div class="detail-table"><div><span>DATE</span><b>${safeText(p.date || p.interval || p.coverage_start)}</b></div><div><span>REGION</span><b>${safeText(p.region)}</b></div><div><span>STATE</span><b>${safeText(p.classification_state)}</b></div><div><span>WITNESSES</span><b>${channels}</b></div><div><span>COORDS</span><b>${Number(p.lat).toFixed(3)}, ${Number(p.lng).toFixed(3)} · ${safeText(p.coordinate_precision)}</b></div><div><span>ROW ID</span><b>${safeText(p.row_id)}</b></div></div>${p.source_url ? `<a class="source-link" target="_blank" rel="noreferrer" href="${p.source_url}">OPEN PRIMARY / SOURCE PAGE ↗</a>` : ''}`;
  }

  function buildEventList(points) {
    const list = $('earth-list');
    list.classList.remove('skeleton-list');
    list.innerHTML = '';
    if (!points.length) { list.innerHTML = '<div class="microcopy">No mappable event rows. Ledger remains available in the repository.</div>'; return; }
    points.forEach(p => {
      const el = document.createElement('div');
      el.className = 'event-card';
      const provisional = (p.source_grade || '').includes('PENDING') || (p.source_grade || '').includes('SERIES_BOUND');
      el.innerHTML = `<div class="row"><b>${safeText(p.site)}</b><small class="badge ${provisional ? 'warn' : ''}">${gradeShort(p.source_grade)}</small></div><div class="row kind"><span>${safeText(p.date || p.row_kind)}</span><small>${safeText(p.region)}</small></div>`;
      el.addEventListener('click', () => { showDetail(p, 'earth'); state.earthSphere?.focus(p); });
      list.appendChild(el);
    });
  }

  function initSpheres() {
    state.earthSphere = new JanusSphere($('earth-globe'), { mode: 'earth', onHover: p => showDetail(p, 'earth'), onSelect: p => showDetail(p, 'earth') }).setPoints(state.earthPoints);
    state.spaceSphere = new JanusSphere($('space-globe'), { mode: 'space', onHover: p => showDetail(p, 'space'), onSelect: p => showDetail(p, 'space') }).setPoints(state.platePoints);
  }

  function setMode(mode) {
    state.mode = mode;
    document.querySelectorAll('.mode').forEach(b => b.classList.toggle('active', b.dataset.mode === mode));
    $('earth-globe').classList.toggle('active', mode === 'earth');
    $('space-globe').classList.toggle('active', mode === 'space');
    $('bridge-view').classList.toggle('active', mode === 'bridge');
    if (mode === 'earth') {
      $('stage-eyebrow').textContent = 'TERRESTRIAL WITNESS LANE'; $('stage-title').textContent = 'Interactive Earth'; $('stage-big').textContent = state.earthPoints.length || '0'; $('stage-small').textContent = 'mapped source-bound observation sites'; $('hud-source').textContent = 'EARTH LEDGER'; $('hud-note').textContent = 'points are evidence locations, not causal links';
    } else if (mode === 'space') {
      $('stage-eyebrow').textContent = 'CELESTIAL SURVEY LANE'; $('stage-title').textContent = 'POSS-I Observation Sphere'; $('stage-big').textContent = state.platePoints.length || '0'; $('stage-small').textContent = 'public plate centers'; $('hud-source').textContent = 'SPACE LEDGER'; $('hud-note').textContent = 'plate centers are survey positions, not anomaly claims';
    } else {
      $('stage-eyebrow').textContent = 'CONTROLLED CROSS-DOMAIN LANE'; $('stage-title').textContent = 'Earth ↔ Space Bridge'; $('stage-big').textContent = 'L2'; $('stage-small').textContent = 'exact-day join gate'; $('hud-source').textContent = 'JOIN CONTRACT'; $('hud-note').textContent = 'join only after compatible time and denominator rights';
    }
  }

  document.querySelectorAll('.mode').forEach(btn => btn.addEventListener('click', () => setMode(btn.dataset.mode)));

  async function load() {
    $('load-state').textContent = 'LOADING';
    const tasks = await Promise.allSettled([fetchJSON(DATA.earth), fetchJSON(DATA.space), fetchJSON(DATA.bridge), fetchText(DATA.plates)]);
    state.earth = tasks[0].status === 'fulfilled' ? tasks[0].value : { rows: fallbackEarthRows };
    state.space = tasks[1].status === 'fulfilled' ? tasks[1].value : null;
    state.bridge = tasks[2].status === 'fulfilled' ? tasks[2].value : null;
    const earthRows = Array.isArray(state.earth?.rows) ? state.earth.rows : fallbackEarthRows;
    state.earthPoints = earthRows.map(rowToPoint).filter(Boolean);
    if (!state.earthPoints.length) state.earthPoints = fallbackEarthRows.map(rowToPoint).filter(Boolean);
    const plateCSV = tasks[3].status === 'fulfilled' ? tasks[3].value : fallbackPlateCSV;
    state.platePoints = parsePlateManifest(plateCSV);
    if (!state.platePoints.length) state.platePoints = parsePlateManifest(fallbackPlateCSV);
    $('earth-count').textContent = `${state.earthPoints.length} MAPPED`;
    buildEventList(state.earthPoints);
    const cohort = state.space?.rows?.find(r => r.row_kind === 'ANALYSIS_COHORT');
    const dataset = state.space?.rows?.find(r => r.row_kind === 'PUBLIC_DATASET_RELEASE');
    $('m-plates').textContent = safeText(dataset?.plate_count || cohort?.plates_in_window);
    $('m-tiles').textContent = safeText(dataset?.tile_rows);
    const candidateCount = dataset?.catalogue_rows || cohort?.candidate_rows_on_in_window_plates;
    $('m-candidates').textContent = candidateCount ? Number(candidateCount).toLocaleString() : '—';
    $('m-days').textContent = safeText(cohort?.unique_observation_dates);
    $('space-count').textContent = `${state.platePoints.length} SKY POSITIONS`;
    if (state.bridge) {
      const l2 = state.bridge.join_levels?.find(x => x.level === 'L2_EXACT_DAY');
      $('bridge-status').textContent = l2?.allowed ? 'L2 EXACT-DAY JOIN: OPEN' : 'L2 EXACT-DAY JOIN: BLOCKED';
      $('bridge-reason').textContent = l2?.current_blocker || l2?.unlock || 'Bridge contract loaded.';
    }
    initSpheres();
    setMode('earth');
    const live = tasks.slice(0, 3).every(t => t.status === 'fulfilled');
    const plateLive = tasks[3].status === 'fulfilled';
    $('load-state').textContent = live && plateLive ? 'LIVE' : 'RESILIENT';
    $('load-state').style.color = live && plateLive ? '#48e59a' : '#ffbd5a';
  }

  load().catch(err => {
    console.error(err);
    $('load-state').textContent = 'DEGRADED'; $('load-state').style.color = '#ff6b6b';
    state.earthPoints = fallbackEarthRows.map(rowToPoint).filter(Boolean);
    state.platePoints = parsePlateManifest(fallbackPlateCSV);
    buildEventList(state.earthPoints); initSpheres(); setMode('earth');
  });
})();