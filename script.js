/* ── NAV ACTIVE ── */
const sections = document.querySelectorAll('section[id], nav');
const navLinks = document.querySelectorAll('nav a');
window.addEventListener('scroll', () => {
  let cur = '';
  document.querySelectorAll('section[id]').forEach(s => {
    if (window.scrollY >= s.offsetTop - 80) cur = s.id;
  });
  navLinks.forEach(a => {
    a.classList.toggle('active', a.getAttribute('href') === '#' + cur);
  });
});

/* ── SLIDER BINDINGS ── */
function bind(id, outId, fmt) {
  const el = document.getElementById(id);
  const out = document.getElementById(outId);
  el.addEventListener('input', () => { out.textContent = fmt(el.value); });
  out.textContent = fmt(el.value);
}
bind('nSectors',   'nVal',  v => v);
bind('badRatio',   'brVal', v => (parseFloat(v)*100).toFixed(0)+'%');
bind('kTargets',   'kVal',  v => v);
bind('clusterProb','cpVal', v => parseFloat(v).toFixed(2));
bind('retryTime',  'rtVal', v => parseFloat(v).toFixed(1)+'s');
bind('startLBA',   'stVal', v => {
  const n = parseInt(document.getElementById('nSectors').value);
  return Math.round(parseInt(v)/100 * n);
});

/* ── SIMULATION ENGINE (FIXED) ── */
/*
 * ROOT CAUSE FIX:
 * Pada graf linear 1D, jarak antara A→B selalu melewati semua sektor di antaranya.
 * Penalti bad sector harus dihitung untuk SETIAP sektor yang DILEWATI, bukan hanya
 * node tujuan. Inilah yang membuat Greedy bisa terjebak: cluster bad sector di
 * antara start dan target A membuat A sebenarnya lebih mahal dari B yang lebih jauh
 * secara spasial tapi lebih bersih jalurnya.
 *
 * GREEDY BUG FIX:
 * Greedy harus melakukan lookup dari posisi SAAT INI (bukan selalu dari index 0).
 * Gunakan distMatrix yang sudah dibangun dari setiap node penting.
 */
function runSimulation() {
  const N  = parseInt(document.getElementById('nSectors').value);
  const br = parseFloat(document.getElementById('badRatio').value);
  const k  = parseInt(document.getElementById('kTargets').value);
  const cp = parseFloat(document.getElementById('clusterProb').value);
  const rt = parseFloat(document.getElementById('retryTime').value);
  const startPct = parseInt(document.getElementById('startLBA').value) / 100;
  const start = Math.round(startPct * (N - 1));

  const SECTOR_MB = 512 / (1024 * 1024);
  const READ_SPEED = 150; const THROUGHPUT = 200;
  const effSpeed = Math.min(READ_SPEED, THROUGHPUT);
  const BASE_COST = SECTOR_MB / effSpeed; // cost per sector (no bad)

  /* ── 1. Generate bad sectors via LOGI BFS ── */
  const rng = seededRng(42);
  const badSet = new Set();
  const status = new Uint8Array(N);
  const bfsQ = [];
  const targetBad = Math.max(1, Math.floor(N * br));

  function plantSeed() {
    for (let tries = 0; tries < 200; tries++) {
      const lba = Math.floor(rng() * N);
      if (!status[lba]) { status[lba]=1; badSet.add(lba); bfsQ.push(lba); return; }
    }
  }
  const numSeeds = Math.max(1, Math.floor(targetBad / 15));
  for (let i = 0; i < numSeeds; i++) plantSeed();

  let qi = 0;
  while (badSet.size < targetBad) {
    if (qi >= bfsQ.length) { plantSeed(); if (qi >= bfsQ.length) break; }
    const lba = bfsQ[qi++];
    for (const nb of [lba-1, lba+1]) {
      if (nb>=0 && nb<N && !status[nb] && rng()<cp) {
        status[nb]=1; badSet.add(nb); bfsQ.push(nb);
        if (badSet.size >= targetBad) break;
      }
    }
  }

  /* ── 2. Pick k targets (good sectors only, well-spread) ── */
  // Bagi N ke dalam k+1 zona, ambil satu sektor good dari tiap zona
  // → menjamin distribusi yang menarik secara visual dan algoritmik
  const targets = [];
  const zoneSize = Math.floor(N / (k + 1));
  const goodPool = [];
  for (let i = 0; i < N; i++) if (!badSet.has(i) && i !== start) goodPool.push(i);
  shuffle(goodPool, rng);

  // Coba ambil satu per zona, fallback ke random
  for (let z = 0; z < k; z++) {
    const lo = (z + 1) * zoneSize - Math.floor(zoneSize * 0.4);
    const hi = (z + 1) * zoneSize + Math.floor(zoneSize * 0.4);
    const zoneGood = goodPool.filter(x => x >= lo && x <= hi && !targets.includes(x));
    if (zoneGood.length) targets.push(zoneGood[Math.floor(rng() * zoneGood.length)]);
  }
  // Isi sisa jika kurang
  for (const x of goodPool) {
    if (targets.length >= k) break;
    if (!targets.includes(x)) targets.push(x);
  }
  targets.sort((a,b) => a-b);

  /* ── 3. Dijkstra dengan akumulasi penalti SEMUA sektor yang dilewati ──
   * w(i → j) = BASE_COST + retry_time jika j adalah bad sector
   * Artinya saat melintasi sebuah jalur dari A ke B, setiap bad sector
   * yang DILEWATI ikut menambah penalti — bukan hanya node tujuan.
   * Ini yang membuat jarak tidak selalu proporsional dengan |A-B|.
   */
  function dijkstra(source) {
    const dist = new Float64Array(N).fill(Infinity);
    dist[source] = 0;
    // Min-heap sederhana: array of [cost, node], diurutkan saat pop
    const pq = [[0.0, source]];
    while (pq.length) {
      // O(n log n) — cukup untuk N≤500
      pq.sort((a,b) => a[0]-b[0]);
      const [cost, u] = pq.shift();
      if (cost > dist[u] + 1e-12) continue;
      for (const nb of [u-1, u+1]) {
        if (nb < 0 || nb >= N) continue;
        // Penalti kena saat MEMASUKI nb (mencakup semua sektor yang dilalui)
        const w = BASE_COST + (badSet.has(nb) ? rt : 0);
        const nc = dist[u] + w;
        if (nc < dist[nb] - 1e-12) { dist[nb] = nc; pq.push([nc, nb]); }
      }
    }
    return dist;
  }

  /* ── 4. Bangun distance matrix antar semua node penting ── */
  const nodes = [start, ...targets]; // index 0 = start, 1..k = targets
  const m = nodes.length;
  const distMatrix = []; // distMatrix[i][j] = waktu dari nodes[i] ke nodes[j]
  for (let i = 0; i < m; i++) {
    const d = dijkstra(nodes[i]);
    distMatrix.push(nodes.map(n => d[n]));
  }

  /* ── 5. Held-Karp DP (exact TSP) ── */
  const FULL = (1 << k) - 1;
  const INF  = Infinity;
  // dp[mask][last] = min cost mengunjungi set(mask), terakhir di target[last]
  const dp  = Array.from({length: 1<<k}, () => new Float64Array(k).fill(INF));
  const par = Array.from({length: 1<<k}, () => new Int16Array(k).fill(-1));

  // Init: dari start (node 0) langsung ke masing-masing target
  for (let i = 0; i < k; i++) dp[1<<i][i] = distMatrix[0][i+1];

  for (let mask = 1; mask <= FULL; mask++) {
    for (let last = 0; last < k; last++) {
      if (!((mask>>last)&1) || dp[mask][last]===INF) continue;
      for (let nxt = 0; nxt < k; nxt++) {
        if ((mask>>nxt)&1) continue;
        const newMask = mask|(1<<nxt);
        // distMatrix[last+1][nxt+1]: dari target[last] ke target[nxt]
        const newCost = dp[mask][last] + distMatrix[last+1][nxt+1];
        if (newCost < dp[newMask][nxt] - 1e-12) {
          dp[newMask][nxt] = newCost;
          par[newMask][nxt] = last;
        }
      }
    }
  }

  // Cari last node terbaik
  let dpTime = INF, bestLast = 0;
  for (let last = 0; last < k; last++) {
    if (dp[FULL][last] < dpTime) { dpTime = dp[FULL][last]; bestLast = last; }
  }

  // Rekonstruksi urutan DP
  const dpOrder = [];
  let rmask = FULL, rcur = bestLast;
  while (rmask) {
    dpOrder.push(rcur);
    const prev = par[rmask][rcur];
    rmask ^= (1<<rcur);
    rcur = prev;
  }
  dpOrder.reverse();
  const dpLBA = dpOrder.map(i => targets[i]);

  /* ── 6. Greedy Nearest Neighbor (FIXED) ──
   * Pada setiap langkah, cari target belum dikunjungi dengan
   * distMatrix[currentNodeIndex][targetIndex] terkecil.
   * currentNodeIndex adalah index di array `nodes` dari posisi saat ini.
   */
  const grVisited = new Set();
  let grTime = 0;
  const grLBA = [];
  let grPosIdx = 0; // index di `nodes`: 0=start, 1..k=targets

  for (let step = 0; step < k; step++) {
    let bestJ = -1, bestD = INF;
    for (let j = 0; j < k; j++) {
      if (grVisited.has(j)) continue;
      const d = distMatrix[grPosIdx][j+1]; // j+1 karena nodes[0]=start
      if (d < bestD - 1e-12) { bestD = d; bestJ = j; }
    }
    grTime += bestD;
    grVisited.add(bestJ);
    grLBA.push(targets[bestJ]);
    grPosIdx = bestJ + 1; // pindah ke target[bestJ] = nodes[bestJ+1]
  }

  /* ── 7. Render canvas ── */
  const canvas = document.getElementById('hddCanvas');
  const ctx = canvas.getContext('2d');
  const DPR = devicePixelRatio || 1;
  canvas.width  = canvas.offsetWidth * DPR;
  canvas.height = 110 * DPR;
  ctx.scale(DPR, DPR);
  const W = canvas.offsetWidth;
  ctx.clearRect(0, 0, W, 110);
  const sw = W / N;

  // ── Sektor background (tengah, y=44-22 height=18) ──
  for (let i = 0; i < N; i++) {
    ctx.fillStyle = badSet.has(i) ? 'rgba(226,75,74,0.5)' : 'rgba(56,90,78,0.18)';
    ctx.fillRect(i*sw, 44, Math.max(sw-0.2, 0.2), 22);
  }

  // ── Divider lines ──
  ctx.strokeStyle = 'rgba(29,158,117,0.12)';
  ctx.lineWidth = 0.5;
  ctx.beginPath(); ctx.moveTo(0,42); ctx.lineTo(W,42); ctx.stroke();
  ctx.beginPath(); ctx.moveTo(0,68); ctx.lineTo(W,68); ctx.stroke();

  // ── Label zona ──
  ctx.font = 'bold 8px monospace';
  ctx.fillStyle = '#9FE1CB'; ctx.fillText('DP', 2, 14);
  ctx.fillStyle = '#F0997B'; ctx.fillText('GR', 2, 108);

  // ── Start marker ──
  ctx.fillStyle = 'rgba(255,152,0,0.9)';
  ctx.fillRect(start*sw, 38, Math.max(sw*2.5, 2), 34);

  // ── Target markers ──
  dpLBA.forEach(lba => {
    ctx.fillStyle = 'rgba(29,158,117,0.85)';
    ctx.fillRect(lba*sw, 38, Math.max(sw*2, 1.5), 34);
  });
  grLBA.forEach((lba, i) => {
    if (grLBA[i] !== dpLBA[i]) {
      ctx.fillStyle = 'rgba(216,90,48,0.6)';
      ctx.fillRect(lba*sw, 38, Math.max(sw*2, 1.5), 34);
    }
  });

  // ── DP path arrows (top zone y=25) ──
  const dpFull = [start, ...dpLBA];
  ctx.strokeStyle = 'rgba(29,158,117,0.9)';
  ctx.lineWidth = 1.8;
  for (let i = 0; i < dpFull.length-1; i++) {
    const x1 = (dpFull[i]+0.5)*sw, x2 = (dpFull[i+1]+0.5)*sw;
    const dir = x2 > x1 ? 1 : -1;
    ctx.beginPath(); ctx.moveTo(x1, 26); ctx.lineTo(x2-dir*5, 26); ctx.stroke();
    ctx.beginPath(); ctx.moveTo(x2-dir*5, 21); ctx.lineTo(x2, 26); ctx.lineTo(x2-dir*5, 31); ctx.stroke();
    ctx.fillStyle = '#9FE1CB'; ctx.font = '8px monospace';
    ctx.fillText(i+1, Math.min(x1,x2)+(Math.abs(x2-x1))/2-3, 18);
  }

  // ── Greedy path arrows (bottom zone y=87) ──
  const grFull = [start, ...grLBA];
  ctx.strokeStyle = 'rgba(216,90,48,0.9)';
  ctx.lineWidth = 1.8;
  for (let i = 0; i < grFull.length-1; i++) {
    const x1 = (grFull[i]+0.5)*sw, x2 = (grFull[i+1]+0.5)*sw;
    const dir = x2 > x1 ? 1 : -1;
    ctx.beginPath(); ctx.moveTo(x1, 84); ctx.lineTo(x2-dir*5, 84); ctx.stroke();
    ctx.beginPath(); ctx.moveTo(x2-dir*5, 79); ctx.lineTo(x2, 84); ctx.lineTo(x2-dir*5, 89); ctx.stroke();
    ctx.fillStyle = '#F0997B'; ctx.font = '8px monospace';
    ctx.fillText(i+1, Math.min(x1,x2)+(Math.abs(x2-x1))/2-3, 103);
  }

  /* ── 8. Tampilkan hasil ── */
  document.getElementById('dpTime').textContent = dpTime.toFixed(4) + 's';
  document.getElementById('grTime').textContent = grTime.toFixed(4) + 's';
  document.getElementById('dpOrder').textContent = '→ LBA: [' + dpLBA.join(' → ') + ']';
  document.getElementById('grOrder').textContent = '→ LBA: [' + grLBA.join(' → ') + ']';
  document.getElementById('dpPath').textContent  = `DP path  : start=${start} → ${dpLBA.join(' → ')}`;
  document.getElementById('grPath').textContent  = `Greedy   : start=${start} → ${grLBA.join(' → ')}`;

  const impr = grTime > 1e-9 ? (grTime - dpTime) / grTime * 100 : 0;
  const badge = document.getElementById('improvBadge');

  // Tampilkan breakdown penalti tiap segmen
  const segBreak = dpLBA.map((lba, i) => {
    const from = i===0 ? start : dpLBA[i-1];
    const lo = Math.min(from, lba)+1, hi = Math.max(from, lba);
    let bad = 0;
    for (let s = lo; s <= hi; s++) if (badSet.has(s)) bad++;
    return `${from}→${lba}(${bad} bad)`;
  }).join(', ');
  document.getElementById('dpPath').textContent += ` | Segmen: [${segBreak}]`;

  if (impr > 0.01) {
    badge.innerHTML = `✅ Hybrid DP lebih cepat <strong>${impr.toFixed(3)}%</strong> dari Greedy &nbsp;|&nbsp; DP: ${dpTime.toFixed(4)}s &nbsp;vs&nbsp; Greedy: ${grTime.toFixed(4)}s`;
    badge.style.background   = 'rgba(29,158,117,0.1)';
    badge.style.borderColor  = 'rgba(29,158,117,0.3)';
    badge.style.color        = 'var(--accent2)';
  } else if (impr < -0.01) {
    badge.innerHTML = `⚠️ Greedy lebih cepat ${Math.abs(impr).toFixed(3)}% — kasus langka (tidak seharusnya terjadi pada DP exact)`;
    badge.style.background   = 'rgba(216,90,48,0.1)';
    badge.style.borderColor  = 'rgba(216,90,48,0.3)';
    badge.style.color        = '#F0997B';
  } else {
    badge.innerHTML = `➖ Kedua algoritma menghasilkan rute berbeda tapi total waktu hampir sama — distribusi bad sector tidak menciptakan "jebakan" pada kondisi ini`;
    badge.style.background   = 'rgba(136,135,128,0.08)';
    badge.style.borderColor  = 'rgba(136,135,128,0.2)';
    badge.style.color        = 'var(--text2)';
  }

  document.getElementById('simOutput').style.display = 'block';
}

/* ── PRESET SCENARIOS ── */
function applyPreset(n, br, k, cp, rt, startPct) {
  document.getElementById('nSectors').value    = n;
  document.getElementById('badRatio').value    = br;
  document.getElementById('kTargets').value    = k;
  document.getElementById('clusterProb').value = cp;
  document.getElementById('retryTime').value   = rt;
  document.getElementById('startLBA').value    = startPct;
  ['nSectors','badRatio','kTargets','clusterProb','retryTime','startLBA'].forEach(id => {
    document.getElementById(id).dispatchEvent(new Event('input'));
  });
  setTimeout(runSimulation, 50);
}

/* ── UTILS ── */
function seededRng(seed) {
  let s = seed;
  return () => {
    s = (s * 1664525 + 1013904223) & 0xffffffff;
    return (s >>> 0) / 4294967296;
  };
}
function shuffle(arr, rng) {
  for (let i = arr.length - 1; i > 0; i--) {
    const j = Math.floor(rng() * (i + 1));
    [arr[i], arr[j]] = [arr[j], arr[i]];
  }
}

/* ── CHARTS ── */
const chartDefaults = {
  responsive: true, maintainAspectRatio: false,
  plugins: { legend: { display: false } },
  elements: { point: { radius: 4 } }
};

/* Chart 1: Improvement per damage level */
new Chart(document.getElementById('chart1'), {
  type: 'bar',
  data: {
    labels: ['5%', '15%', '30%'],
    datasets: [{
      label: 'Avg Improvement (%)',
      data: [2.1, 6.8, 12.07],
      backgroundColor: ['rgba(29,158,117,0.6)', 'rgba(29,158,117,0.75)', 'rgba(29,158,117,0.9)'],
      borderColor: ['#1D9E75', '#1D9E75', '#1D9E75'],
      borderWidth: 1,
      borderRadius: 4
    }]
  },
  options: {
    ...chartDefaults,
    scales: {
      x: { ticks: { color: '#8fb09a', font: { size: 12 } }, grid: { color: 'rgba(29,158,117,0.08)' } },
      y: {
        ticks: { color: '#8fb09a', font: { size: 11 }, callback: v => v + '%' },
        grid: { color: 'rgba(29,158,117,0.08)' },
        min: 0, max: 14
      }
    },
    plugins: {
      ...chartDefaults.plugins,
      tooltip: { callbacks: { label: ctx => ctx.raw.toFixed(2) + '%' } }
    }
  }
});

/* Chart 2: CPU time vs k */
new Chart(document.getElementById('chart2'), {
  type: 'line',
  data: {
    labels: [3, 5, 8, 10, 12, 15],
    datasets: [
      {
        label: 'Hybrid DP',
        data: [0.8, 2.1, 18.4, 89.3, 210.5, 487.2],
        borderColor: '#1D9E75', backgroundColor: 'rgba(29,158,117,0.1)',
        borderWidth: 2, fill: true,
        borderDash: []
      },
      {
        label: 'Greedy',
        data: [0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
        borderColor: '#D85A30', backgroundColor: 'rgba(216,90,48,0.08)',
        borderWidth: 2, fill: true,
        borderDash: [4, 4]
      }
    ]
  },
  options: {
    ...chartDefaults,
    scales: {
      x: { ticks: { color: '#8fb09a', font: { size: 12 } }, grid: { color: 'rgba(29,158,117,0.08)' } },
      y: {
        ticks: { color: '#8fb09a', font: { size: 11 }, callback: v => v + 'ms' },
        grid: { color: 'rgba(29,158,117,0.08)' }, min: 0
      }
    },
    plugins: {
      ...chartDefaults.plugins,
      tooltip: { callbacks: { label: ctx => ctx.dataset.label + ': ' + ctx.raw.toFixed(1) + 'ms' } }
    }
  }
});

/* Chart 3: Recovery time vs bad ratio */
const brLabels = ['5%', '10%', '15%', '20%', '25%', '30%'];
new Chart(document.getElementById('chart3'), {
  type: 'line',
  data: {
    labels: brLabels,
    datasets: [
      {
        label: 'Hybrid DP',
        data: [385.2, 421.8, 458.4, 510.6, 562.3, 620.1],
        borderColor: '#1D9E75', backgroundColor: 'rgba(29,158,117,0.1)',
        borderWidth: 2, fill: false
      },
      {
        label: 'Greedy',
        data: [392.1, 435.7, 490.2, 553.1, 618.4, 703.8],
        borderColor: '#D85A30', backgroundColor: 'rgba(216,90,48,0.08)',
        borderWidth: 2, fill: false, borderDash: [4, 4]
      }
    ]
  },
  options: {
    ...chartDefaults,
    scales: {
      x: { ticks: { color: '#8fb09a', font: { size: 12 } }, grid: { color: 'rgba(29,158,117,0.08)' } },
      y: {
        ticks: { color: '#8fb09a', font: { size: 11 }, callback: v => v + 's' },
        grid: { color: 'rgba(29,158,117,0.08)' }
      }
    }
  }
});

/* Chart 4: Improvement % vs bad ratio */
new Chart(document.getElementById('chart4'), {
  type: 'line',
  data: {
    labels: brLabels,
    datasets: [{
      label: 'Improvement (%)',
      data: [1.8, 3.2, 6.5, 7.7, 9.1, 12.07],
      borderColor: '#7F77DD', backgroundColor: 'rgba(127,119,221,0.12)',
      borderWidth: 2.5, fill: true, tension: 0.3,
      pointBackgroundColor: '#7F77DD', pointRadius: 5
    }]
  },
  options: {
    ...chartDefaults,
    scales: {
      x: { ticks: { color: '#8fb09a', font: { size: 12 } }, grid: { color: 'rgba(29,158,117,0.08)' } },
      y: {
        ticks: { color: '#8fb09a', font: { size: 11 }, callback: v => v + '%' },
        grid: { color: 'rgba(29,158,117,0.08)' }, min: 0, max: 14
      }
    },
    plugins: {
      ...chartDefaults.plugins,
      tooltip: { callbacks: { label: ctx => 'Improvement: ' + ctx.raw.toFixed(2) + '%' } }
    }
  }
});
