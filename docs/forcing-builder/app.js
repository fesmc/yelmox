// Forcing builder: a view on one run's configuration (a base par file plus
// runme -p overrides, resolved against the forcing library and the schema
// defaults). Left: the inputs as cards (time, series, terms, snapshots,
// transforms, sea level and orbit). Right: the equation chain of the forcing
// evaluated at a chosen time, the series over the run, and the run written out
// (the par file, the equation as text, the -p list and the runme command).
// Presets: the cases of the docs figures. Data: forcing.json
// (scripts/forcing-builder/export.jl). The series formulas mirror fesm-utils
// tsgen (TsGen.jl); the equations follow docs/_forcing-terms.md section 2.

import * as Plot from "https://cdn.jsdelivr.net/npm/@observablehq/plot@0.6.16/+esm";
import katex from "https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/katex.mjs";

const root = document.getElementById("fb");
let D = null;                       // forcing.json
let TEMPLATES = [];                 // [regex, template] from forcing.json (series_, term_, snapshot_)
// state: the base, the overrides "g.k" -> value, the -p tokens outside the
// exported groups (passed through), the open details, the preset loaded (and
// its override key, to tell when it was modified), the evaluation time, the tab
const S = { base: null, ov: new Map(), extra: [], open: new Set(), runopts: "-r", outdir: "", preset: null, presetKey: "", tsel: null, tab: "par" };

const FIELDS = ["tas", "pr", "to", "so"];
const FSYM = { tas: "T^{\\mathrm{sl}}", pr: "P^{\\mathrm{sl}}", to: "\\Theta", so: "S" };
const GAIN = { tas: "f_ta", to: "f_to", so: "f_so" };
const GREEK = { alpha: "\\alpha", beta: "\\beta", gamma: "\\gamma", delta: "\\delta", epsilon: "\\epsilon", zeta: "\\zeta", eta: "\\eta", theta: "\\theta", kappa: "\\kappa", lambda: "\\lambda", mu: "\\mu", nu: "\\nu", xi: "\\xi", rho: "\\rho", sigma: "\\sigma", tau: "\\tau", phi: "\\phi", chi: "\\chi", psi: "\\psi", omega: "\\omega" };
const sub = (n) => "\\mathrm{" + String(n).replace(/_/g, "\\_") + "}";
const sym = (n) => GREEK[n] || sub(n);

// ---------------------------------------------------------------- resolution

function templateOf(g) {
  if (D.defaults[g]) return g;
  for (const [re, t] of TEMPLATES) if (re.test(g)) return t;
  return null;
}
const base = () => D.bases[S.base];
function groupSource(g) {
  if (base().groups[g]) return "base";
  if (D.library.groups[g]) return "library";
  if (templateOf(g)) return "new";
  return null;
}
// a named group that exists (in the par file or the library), not only as a template
const groupExists = (g) => ["base", "library"].includes(groupSource(g));
function groupDef(g) {
  const src = groupSource(g);
  if (src === "base") return base().groups[g];
  if (src === "library") return D.library.groups[g];
  return null;
}
// The schema of a group: its keys in order, with comment and the reference value.
function schemaKeys(g) {
  const t = templateOf(g), def = groupDef(g);
  const order = [];
  if (t) order.push(...D.defaults[t].order);
  if (def) for (const k of def.order) if (!order.includes(k)) order.push(k);
  return order;
}
function keyInfo(g, k) {
  const t = templateOf(g), def = groupDef(g), src = groupSource(g);
  const tk = t ? D.defaults[t].keys[k] : null;
  const dk = def ? def.keys[k] : null;
  if (!tk && !dk) return null;
  const ref = dk ? dk.v : tk.v;
  const refsrc = dk ? src : "default";
  const comment = (dk && dk.c) ? dk.c : (tk ? tk.c : "");
  const pre = (dk && dk.pre) ? dk.pre : (tk ? tk.pre : "");
  const id = g + "." + k;
  if (S.ov.has(id)) return { v: S.ov.get(id), src: "override", ref, refsrc, comment, pre };
  return { v: ref, src: refsrc, ref, refsrc, comment, pre };
}
const eff = (g, k, d = undefined) => { const i = keyInfo(g, k); return i ? i.v : d; };
const asList = (v) => Array.isArray(v) ? v.map(String).filter(s => s !== "") : (v === "" || v == null ? [] : [String(v)]);
const same = (a, b) => JSON.stringify(a) === JSON.stringify(b);
const sameValue = (a, b) => (Array.isArray(a) || Array.isArray(b)) ? same(asList(a), asList(b)) : same(a, b);

function setOv(g, k, v) {
  const i = keyInfo(g, k);
  const id = g + "." + k;
  if (i && sameValue(v, i.ref)) S.ov.delete(id); else S.ov.set(id, v);
  render();
}

// the timeline group ([ctrl] run_step)
const tgroup = () => { const r = eff("ctrl", "run_step", "ctrl"); return groupSource(r) ? r : "ctrl"; };
const climate = () => eff("comps", "climate", "snapshots");
// the named groups of a kind in the par file and the library (prefix stripped)
const named = (prefix) => [...new Set([...Object.keys(base().groups), ...Object.keys(D.library.groups)])].filter(g => g.startsWith(prefix)).map(g => g.slice(prefix.length)).sort();

// ------------------------------------------------------------------- helpers

function h(tag, attrs = {}, ...children) {
  const el = document.createElement(tag);
  for (const [a, v] of Object.entries(attrs)) {
    if (a === "class") el.className = v;
    else if (a === "html") el.innerHTML = v;
    else if (a.startsWith("on")) el.addEventListener(a.slice(2), v);
    else if (v !== null && v !== undefined) el.setAttribute(a, v);
  }
  for (const c of children.flat()) if (c != null) el.append(c.nodeType ? c : document.createTextNode(String(c)));
  return el;
}
const fmtNum = (x) => {
  if (typeof x !== "number") return String(x);
  if (Number.isInteger(x) && Math.abs(x) >= 1000 && x % 1000 === 0) return (x / 1000) + "e3";
  return Number.isInteger(x) ? String(x) : String(+x.toPrecision(6));
};
const fmtK = (x) => fmtNum(x).replace(/e3$/, "\\,\\mathrm{k}");
function tex(s, display = true) {
  return h("div", { class: display ? "fb-eq" : "fb-eq-inline", html: katex.renderToString(s, { displayMode: display, throwOnError: false }) });
}
const num = (s) => { const x = Number(String(s).replace(/d/i, "e")); return Number.isFinite(x) ? x : NaN; };

// Known option lists, by "group.key" (or "template.key"); the series method
// list comes from its comment (const|ramp-slope|...|derived).
const OPTIONS = {
  "ctrl.tstep_method": ["const", "rel", "cal"], "spinup.tstep_method": ["const", "rel", "cal"], "transient.tstep_method": ["const", "rel", "cal"],
  "comps.surface_method": ["chion", "smbpal", "smb_simple", "climate"],
  "sim.init_method": ["none", "equil", "recon", "recon_ref"],
  "barysealevel.method": ["const", "file", "fastiso"],
  "orbit.method": ["forcing_time", "fixed_time", "fixed_orbit"], "ice_bnd.method": ["file", "scaled"],
  "snapshots.combine": ["anomaly", "absolute"], "snapshot.tas_season": ["", "jja", "djf"],
  "series.from": ["", "air_anomaly", "obm_north", "obm_south"],
  "term.pattern": ["snapshots", "uniform", "replace"],
  "taper.outside": ["ref", "mean"],
  "tm_1D.method": ["const", "file", "times"], "tm_2D.method": ["const", "file", "times"], "tm_2Dsm.method": ["const", "file", "times"], "tm_rst.method": ["const", "file", "times"],
};
function optionsOf(g, k, info) {
  const t = templateOf(g);
  if (OPTIONS[g + "." + k]) return OPTIONS[g + "." + k];
  if (t && OPTIONS[t + "." + k]) return OPTIONS[t + "." + k];
  const m = (info.comment || "").split("\n")[0].match(/^([\w.-]+(\|[\w.-]+)+)/);
  if (m) return m[1].split("|");
  return null;
}

// A labelled widget of one key: checkbox, select, number or text, the source
// badge and the comment as help. opts: {label, options, labels: {value: label}}.
function field(g, k, opts = {}) {
  const info = keyInfo(g, k);
  if (!info) return h("div", { class: "fb-field fb-missing" }, `${g}.${k}: not in the schema`);
  const v = info.v, id = g + "." + k;
  let input;
  const options = opts.options || optionsOf(g, k, info);
  const labels = opts.labels || {};
  if (typeof v === "boolean") {
    input = h("input", { type: "checkbox", onchange: (e) => setOv(g, k, e.target.checked) });
    input.checked = v;
  } else if (options && !Array.isArray(v)) {
    input = h("select", { onchange: (e) => setOv(g, k, e.target.value) },
      ...options.map(o => h("option", { value: o }, labels[o] ?? o)), ...(options.includes(String(v)) ? [] : [h("option", { value: v }, String(v))]));
    input.value = String(v);
  } else if (typeof v === "number") {
    input = h("input", { type: "text", class: "fb-num", value: fmtNum(v), onchange: (e) => { const x = num(e.target.value); if (Number.isFinite(x)) setOv(g, k, x); else render(); } });
  } else if (Array.isArray(v)) {
    const isnum = v.every(x => typeof x === "number");
    input = h("input", { type: "text", class: "fb-text", value: v.map(fmtNum).join(", "),
      onchange: (e) => { const items = e.target.value.split(/[,\s]+/).filter(s => s !== ""); setOv(g, k, isnum ? items.map(num) : items); } });
  } else {
    input = h("input", { type: "text", class: "fb-text", value: String(v), onchange: (e) => setOv(g, k, e.target.value) });
  }
  return fieldRow(g, k, info, input, opts.label);
}
function fieldRow(g, k, info, input, label) {
  const badge = h("span", { class: "fb-src fb-src-" + info.src, title: `value from the ${info.src}` + (info.src === "override" ? ` (${info.refsrc}: ${JSON.stringify(info.ref)})` : "") }, info.src);
  const help = (info.comment || "").split("\n")[0];
  return h("div", { class: "fb-field" },
    h("label", { title: info.comment, for: g + "." + k }, h("code", {}, label || k)),
    input, badge, help ? h("span", { class: "fb-help" }, help) : null);
}
function fields(g, keys, opts = {}) { return h("div", { class: "fb-fields" }, ...keys.map(k => field(g, k, opts[k] || {}))); }

// A vector-of-strings key as checkboxes over a fixed set (the fields of a term).
function checkList(g, k, options, marks = {}) {
  const info = keyInfo(g, k);
  if (!info) return h("div", { class: "fb-field fb-missing" }, `${g}.${k}: not in the schema`);
  const cur = asList(info.v);
  const input = h("span", { class: "fb-checks" }, ...options.map(o => {
    const cb = h("input", { type: "checkbox", onchange: (e) => setOv(g, k, options.filter(x => x === o ? e.target.checked : cur.includes(x))) });
    cb.checked = cur.includes(o);
    return h("label", { title: marks[o] || "" }, cb, o + (marks[o] ? " ✗" : ""));
  }));
  return fieldRow(g, k, info, input);
}

// A list of names (a vector-of-strings key): chips with remove, and an add control.
function namesEditor(g, k, candidates, label) {
  const names = asList(eff(g, k, ""));
  const cands = candidates.filter(c => !names.includes(c));
  const sel = h("select", {}, h("option", { value: "" }, "add…"), ...cands.map(c => h("option", { value: c }, c)));
  const txt = h("input", { type: "text", class: "fb-num", placeholder: "new name" });
  const add = (n) => { if (n && !names.includes(n)) setOv(g, k, [...names, n]); };
  sel.addEventListener("change", () => add(sel.value));
  txt.addEventListener("change", () => add(txt.value.trim()));
  return h("div", { class: "fb-names" },
    h("code", {}, label || `${g}.${k}`), " = ",
    ...names.map(n => h("span", { class: "fb-chip" }, n, h("button", { class: "fb-x", title: "remove", onclick: () => setOv(g, k, names.filter(x => x !== n)) }, "×"))),
    sel, txt);
}

// ------------------------------------------------------------------ timeline

function timeline() {
  const T = tgroup();
  const method = String(eff(T, "tstep_method", "const"));
  const t0 = num(eff(T, "time_init", 0)), t1 = num(eff(T, "time_end", 0));
  const tc = num(eff(T, "tstep_const", 0));
  const cal = eff("ctrl", "calendar", false), cref = num(eff("ctrl", "calendar_ref", 1950));
  const n = 400;
  const t = Array.from({ length: n + 1 }, (_, i) => t0 + (t1 - t0) * i / n);
  const tf = t.map(x => method === "const" ? tc : (method === "cal" ? x - cref : x));
  return { T, method, t0, t1, tc, cal, cref, t, tf, tf0: tf[0], tfmin: Math.min(...tf), tfmax: Math.max(...tf) };
}

// ----------------------------------------------------------------- tsgen (JS)

const FEEDBACK = ["exp", "PI42", "H312b", "H312PID", "H321PID", "PID1"];
function seriesOf(file, varname) {
  const s = D.series[file];
  if (!s) return null;
  if (s.kind === "netcdf") { const v = s.vars[varname]; return v ? { time: v.time, f: v.f, header: `${file}:${varname}`, n: v.n } : null; }
  return s;
}
function interp(ser, t) {
  const T = ser.time, F = ser.f, n = T.length;
  if (!n) return NaN;
  if (t <= T[0]) return F[0]; if (t >= T[n - 1]) return F[n - 1];
  let lo = 0, hi = n - 1;
  while (hi - lo > 1) { const m = (lo + hi) >> 1; if (T[m] <= t) lo = m; else hi = m; }
  const w = (t - T[lo]) / (T[hi] - T[lo]);
  return F[lo] * (1 - w) + F[hi] * w;
}
function seriesMean(ser, tw) {
  if (tw.length === 1) return interp(ser, tw[0]);
  const lo = Math.min(...tw), hi = Math.max(...tw);
  const sel = ser.f.filter((f, i) => ser.time[i] >= lo && ser.time[i] <= hi);
  return sel.length ? sel.reduce((a, b) => a + b, 0) / sel.length : interp(ser, (lo + hi) / 2);
}
function seriesPar(g) {
  const p = {};
  for (const k of schemaKeys(g)) p[k] = eff(g, k);
  // normalisation is active only when both norm_t0 and norm_t1 are set in the group itself
  const def = groupDef(g);
  p.norm = !!(def && def.keys.norm_t0 && def.keys.norm_t1) || S.ov.has(g + ".norm_t0") || S.ov.has(g + ".norm_t1");
  return p;
}
function seriesFile(p) {
  if (p.method !== "series" || !p.series_file) return null;
  let ser = seriesOf(p.series_file, p.series_var);
  if (!ser) return null;
  if (p.norm) {
    const v0 = seriesMean(ser, asList(p.norm_t0).map(num)), v1 = seriesMean(ser, asList(p.norm_t1).map(num));
    ser = { ...ser, f: ser.f.map(f => (f - v0) / (v1 - v0)), norm: [v0, v1] };
  }
  return ser;
}
function evalOpen(p, te) {
  const fmin = num(p.f_min), fmax = num(p.f_max), s = num(p.df_sign) >= 0 ? 1 : -1;
  const tau = Math.max(0, te - num(p.dt_init));
  const fstart = s > 0 ? fmin : fmax;
  const clamp = (x) => Math.min(Math.max(x, fmin), fmax);
  switch (p.method) {
    case "const": case "ramp-slope": return clamp(fstart + s * num(p.df_dt_max) * tau);
    case "ramp-time": return clamp(fstart + s * Math.abs(fmax - fmin) / num(p.dt_ramp) * tau);
    case "ramp-time-step": {
      const dr = num(p.dt_ramp), dc = num(p.dt_conv), fc = num(p.f_conv), fext = s > 0 ? fmax : fmin;
      if (tau <= dr) return fstart + (fext - fstart) * (tau / dr);
      if (tau <= dr + dc) return fext + (fc - fext) * ((tau - dr) / dc);
      return fc;
    }
    case "sin": { const T = te <= num(p.dt_init) ? 0 : te; return 0.5 * (fmax + fmin) + 0.5 * (fmax - fmin) * Math.sin(2 * Math.PI * T / num(p.dt_ramp)); }
    default: return fstart;
  }
}
function applyWindow(p, t, f) {
  const t0 = num(p.win_t0), t1 = num(p.win_t1);
  if (t1 <= t0 || t <= t0 || t >= t1) return f;
  const ref = num(p.win_ref);
  if (p.win_above && f <= ref) return f;
  const taper = num(p.win_taper);
  const w = taper > 0 ? Math.min(1, (t - t0) / taper, (t1 - t) / taper) : 1;
  return f + w * (num(p.win_fac) - 1) * (f - ref);
}
// The series [series_<name>] on the run's timeline (tl): {p, ser, f}; a derived
// series is set by the model (f = null); feedback methods: the envelope.
function seriesOnRun(name, tl) {
  const g = "series_" + name;
  if (!groupSource(g)) return { p: null, ser: null, f: null, derived: false, missing: true };
  const p = seriesPar(g);
  if (p.method === "derived") return { p, ser: null, f: null, derived: true, missing: false };
  const ser = seriesFile(p);
  const f = tl.tf.map(tf => applyWindow(p, tf, p.method === "series" ? (ser ? interp(ser, tf) : NaN) : evalOpen(p, tf - tl.tf0)));
  return { p, ser, f, derived: false, missing: false };
}
// The weight of a term over the run: ones for series = "", the series values, or null (derived, missing).
function termWeight(t, tl, snames) {
  if (t.series === "") return tl.t.map(() => 1);
  if (!snames.includes(t.series)) return null;
  const s = seriesOnRun(t.series, tl);
  return s.f && s.f.some(Number.isFinite) ? s.f : null;
}

// ------------------------------------------------------------------- figures

// Plot widths: the right column (the run) and the cards of the left column;
// the SVGs scale down with their column (app.css).
const WR = 700, HR = 200, WL = 560, HL = 120;
function linePlot(series, { x = "t [yr]", y = "", rules = [], title = "", vrules = [], band = null, width = WR, height = HR, tmark = null }) {
  // series: [{label, t:[], f:[]}]
  const data = series.flatMap(s => s.t.map((t, i) => ({ t, f: s.f[i], label: s.label })));
  const marks = [Plot.ruleY([0], { stroke: "#999" })];
  if (band) marks.push(Plot.rectX([band], { x1: "x0", x2: "x1", fill: "#ffd", stroke: "none" }));
  marks.push(Plot.line(data, { x: "t", y: "f", stroke: "label", strokeWidth: 1.5 }));
  for (const r of rules) marks.push(Plot.ruleY([r.y], { stroke: r.stroke || "#c33", strokeDasharray: "4 3" }), Plot.text([{ t: data.length ? data[0].t : 0, f: r.y }], { x: "t", y: "f", text: () => r.label, dy: -6, textAnchor: "start", fill: r.stroke || "#c33" }));
  for (const r of vrules) marks.push(Plot.ruleX([r.x], { stroke: r.stroke || "#36c", strokeDasharray: "4 3" }), Plot.text([{ t: r.x }], { x: "t", text: () => r.label, frameAnchor: "top", dx: 4, textAnchor: "start", fill: r.stroke || "#36c" }));
  if (tmark !== null) marks.push(Plot.ruleX([tmark], { stroke: "#555", strokeWidth: 1 }));
  const ys = data.map(d => d.f).filter(Number.isFinite).concat(rules.map(r => r.y));
  let domain;
  if (ys.length) { const lo = Math.min(...ys), hi = Math.max(...ys); if (hi - lo < 1e-9) domain = [lo - 1, hi + 1]; }
  const el = Plot.plot({ width, height, marginLeft: 55, title: title || undefined, x: { label: x }, y: { label: y, grid: true, domain },
    color: series.length > 1 ? { legend: true } : undefined, marks });
  return h("figure", { class: "fb-fig" }, el);
}
const smallPlot = (series, opts = {}) => linePlot(series, { width: WL, height: HL, ...opts });

// ------------------------------------------------------------------- layout

// Left: a section of inputs (uppercase head with a hint). Right: a numbered
// block of the equation chain or of the run.
function section(id, title, hint, ...content) {
  return h("section", { class: "fb-section", id: "fb-" + id }, h("h2", {}, title, hint ? h("span", { class: "fb-hint" }, hint) : null), ...content);
}
function eqblock(title, ...content) { return h("div", { class: "fb-eqblock" }, h("h3", {}, title), ...content); }
function details(id, title, ...content) {
  const d = h("details", { class: "fb-details" }, h("summary", {}, title), ...content);
  if (S.open.has(id)) d.open = true;
  d.addEventListener("toggle", () => { if (d.open) S.open.add(id); else S.open.delete(id); });
  return d;
}
function note(text, cls = "fb-note") { return h("p", { class: cls }, text); }
function warn(text) { return note(text, "fb-warn"); }
function vals(text) { return h("div", { class: "fb-vals" }, text); }
function card(head, src, ...content) {
  return h("div", { class: "fb-card" }, h("div", { class: "fb-cardhead" }, h("code", {}, head), src ? h("span", { class: "fb-src fb-src-" + src }, src) : null), ...content);
}

// The model time at which the equations are evaluated (the slider), within the run.
function evalTime(tl) {
  const lo = Math.min(tl.t0, tl.t1), hi = Math.max(tl.t0, tl.t1);
  return S.tsel === null || !Number.isFinite(S.tsel) ? tl.t0 : Math.min(Math.max(S.tsel, lo), hi);
}
const tauAt = (tl, t) => tl.method === "const" ? tl.tc : (tl.method === "cal" ? t - tl.cref : t);
const atT = (tl, f, t) => f ? interp({ time: tl.t, f }, t) : NaN;
const fmtV = (x) => Number.isFinite(x) ? fmtNum(+x.toPrecision(4)) : "–";

// ------------------------------------------------------------------ top bar

function topBar() {
  const sel = h("select", { onchange: (e) => { S.base = e.target.value; S.ov.clear(); S.extra = []; S.preset = null; S.tsel = null; render(); } },
    ...Object.keys(D.bases).map(b => h("option", { value: b }, b)));
  sel.value = S.base;
  const dom = eff("domain", "name", "?"), grid = eff("domain", "grid_hub", eff("domain", "grid_ice", ""));
  const state = presetState();
  const presets = (D.presets || []).map(p => h("button", { class: "fb-preset" + (S.preset === p.name ? (state === "on" ? " on" : " mod") : ""), title: `${D.bases[p.base].file} -p ${p.p.join(" ")}`, onclick: () => loadPreset(p) }, p.name));
  return h("div", { class: "fb-top" },
    h("div", { class: "fb-toprow" },
      h("label", {}, h("code", {}, "-n"), " base par file ", sel),
      h("span", { class: "fb-help" }, `${base().file}: domain ${dom}${grid ? ", grid " + grid : ""}, climate backend ${climate()}`),
      h("button", { class: "fb-btn", onclick: () => { S.ov.clear(); S.extra = []; S.preset = null; render(); } }, "reset to base")),
    h("div", { class: "fb-toprow" }, h("span", { class: "fb-hint" }, "presets: the cases of the docs figures"), h("div", { class: "fb-presets" }, ...presets),
      S.preset ? h("span", { class: "fb-help" }, state === "on" ? `${S.preset}` : `${S.preset}, modified`) : null));
}
function ovKey() { return JSON.stringify([[...S.ov.entries()].sort(), S.extra]); }
function presetState() { return S.preset && S.presetKey === ovKey() ? "on" : "mod"; }
function loadPreset(p) {
  S.base = p.base; S.ov.clear(); S.extra = []; S.tsel = null;
  parseP(p.p.join(" "));
  S.preset = p.name; S.presetKey = ovKey();
  render();
}

// ------------------------------------------------------- left: run control

function secRun(tl) {
  const T = tl.T;
  const stepOpts = ["ctrl", ...["spinup", "transient"].filter(g => groupSource(g) === "base")];
  const content = [
    stepOpts.length > 1 ? fields("ctrl", ["run_step"], { run_step: { options: stepOpts } }) : null,
    fields(T, ["tstep_method", "tstep_const", "time_init", "time_end", "dtt"]),
    fields("ctrl", ["calendar", "calendar_ref", "restart"]),
    note("tstep_method sets the forcing time τ (equation 0): const holds the forcing at tstep_const while the model time advances; rel follows the model time (paleo records); cal: calendar years against calendar_ref."),
  ];
  if (keyInfo(T, "time_ref")) content.push(details("esm-periods", "ESM periods", fields(T, ["time_ref", "time_esm_ref", "time_hist", "time_proj", "clim_var", "clim_seed"])));
  content.push(details("components", "Components, output, initial state",
    fields("comps", ["climate", "surface_method", "dt_clim", "with_ice_sheet", "with_isostasy", "with_climate", "with_surface", "with_shelf"], { climate: { options: [String(eff("comps", "climate"))] } }),
    note("The climate backend is fixed by the base: pick another base for another backend."),
    h("div", { class: "fb-fields" }, field("tm_1D", "dt", { label: "tm_1D.dt" }), field("tm_2D", "dt", { label: "tm_2D.dt" }), field("tm_rst", "dt", { label: "tm_rst.dt" })),
    fields("sim", ["init_method", "relax", "opt"]),
    ...(() => { const im = String(eff("sim", "init_method", "none")); return im.startsWith("recon") ? [fields("sim", ["recon_path", "recon_var", "recon_time", "recon_regions"])] : im === "equil" ? [fields("sim", ["init_equil_time"])] : []; })()));
  return section("run", "Time", `[${T}] defines the forcing time τ`, h("div", { class: "fb-card" }, ...content));
}

// ------------------------------------------------------------ left: series

const METHOD_EQ = {
  "const": (p, a) => `${a} = \\min\\!\\big(\\max(f_{\\min} + s\\,\\dot f_{\\max}\\,\\tau_e,\\ f_{\\min}),\\ f_{\\max}\\big)`,
  "ramp-slope": (p, a) => `${a} = \\min\\!\\big(\\max(f_{\\mathrm{start}} + s\\,\\dot f_{\\max}\\,\\tau_e,\\ f_{\\min}),\\ f_{\\max}\\big)`,
  "ramp-time": (p, a) => `${a} = \\min\\!\\big(\\max(f_{\\mathrm{start}} + s\\,\\tfrac{|f_{\\max}-f_{\\min}|}{\\Delta t_{\\mathrm{ramp}}}\\,\\tau_e,\\ f_{\\min}),\\ f_{\\max}\\big)`,
  "ramp-time-step": (p, a) => `${a} = \\begin{cases} f_{\\mathrm{start}} + (f_{\\mathrm{ext}}-f_{\\mathrm{start}})\\,\\tau_e/\\Delta t_{\\mathrm{ramp}} & \\tau_e \\le \\Delta t_{\\mathrm{ramp}} \\\\ f_{\\mathrm{ext}} + (f_{\\mathrm{conv}}-f_{\\mathrm{ext}})\\,(\\tau_e-\\Delta t_{\\mathrm{ramp}})/\\Delta t_{\\mathrm{conv}} & \\tau_e \\le \\Delta t_{\\mathrm{ramp}}+\\Delta t_{\\mathrm{conv}} \\\\ f_{\\mathrm{conv}} & \\text{after} \\end{cases}`,
  "sin": (p, a) => `${a} = \\tfrac{f_{\\max}+f_{\\min}}{2} + \\tfrac{f_{\\max}-f_{\\min}}{2}\\,\\sin\\!\\big(2\\pi\\, t_e/\\Delta t_{\\mathrm{ramp}}\\big)`,
  "series": (p, a) => (p.norm ? `${a} = \\frac{F(\\tau) - \\bar F[t_0]}{\\bar F[t_1] - \\bar F[t_0]}` : `${a} = F(\\tau)`) + `,\\quad F = \\texttt{${String(p.series_file).replace(/^input\//, "").replace(/_/g, "\\_")}}`,
};
const METHOD_KEYS = {
  "const": ["f_min", "f_max", "df_dt_max", "df_sign", "dt_init"],
  "ramp-slope": ["f_min", "f_max", "df_dt_max", "df_sign", "dt_init"],
  "ramp-time": ["f_min", "f_max", "dt_ramp", "df_sign", "dt_init"],
  "ramp-time-step": ["f_min", "f_max", "dt_ramp", "dt_conv", "f_conv", "df_sign", "dt_init"],
  "sin": ["f_min", "f_max", "dt_ramp", "dt_init"],
  "series": ["series_file", "series_var", "series_time_name", "norm_t0", "norm_t1"],
};
const FEEDBACK_KEYS = ["f_min", "f_max", "df_dt_max", "df_sign", "dt_init", "dt_ave", "tol", "with_kill"];
const DERIVED_KEYS = { air_anomaly: ["factor", "region"], obm_north: ["anomaly"], obm_south: ["anomaly"] };

// The definition of a derived series (set by the model each step).
function derivedEq(name, p) {
  const a = sym(name) + "(\\tau)";
  const from = String(p.from || "");
  if (from === "air_anomaly") return `${a} = ${fmtNum(num(p.factor))}\\,\\overline{\\big(T^{\\mathrm{sl}}_{\\mathrm{ann}} - T^{\\mathrm{sl}}_{\\mathrm{ann},0}\\big)}^{\\,${p.region === "all" ? "x" : sub(p.region)}}`;
  if (from === "obm_north" || from === "obm_south") { const b = from === "obm_north" ? "n" : "s"; return `${a} = \\theta_{${b}}` + (p.anomaly ? ` - \\theta_{${b},\\mathrm{init}}` : ""); }
  return `${a} = \\text{?}\\qquad(\\texttt{from} = ${sub(from || "\"\"")})`;
}
// The equation of a series (block 1 on the right).
function seriesEq(name, p) {
  const a = sym(name) + "(\\tau)", method = String(p.method);
  if (method === "derived") return derivedEq(name, p) + `\\qquad\\text{(derived: set by the model each step)}`;
  if (FEEDBACK.includes(method)) return `\\dot ${sym(name)} = \\text{${method}}(\\dot V_{\\mathrm{ice}}),\\quad |\\dot ${sym(name)}| \\le ${fmtNum(p.df_dt_max)},\\quad ${sym(name)} \\in [${fmtNum(p.f_min)},\\ ${fmtNum(p.f_max)}]`;
  return METHOD_EQ[method] ? METHOD_EQ[method](p, a) : `${a} = \\text{?}`;
}

function seriesCard(name, tl) {
  const g = "series_" + name;
  const src = groupSource(g);
  if (!src) return h("div", { class: "fb-card" }, warn(`[${g}]: no such group and no series template`));
  const { p, ser, f, derived } = seriesOnRun(name, tl);
  const method = String(p.method);
  const head = [h("div", { class: "fb-cardhead" }, h("span", { class: "fb-sym", html: katex.renderToString(sym(name) + "(\\tau)", { throwOnError: false }) }), h("code", {}, `[${g}]`), h("span", { class: "fb-src fb-src-" + src }, src === "new" ? "new (template series)" : src),
    h("button", { class: "fb-x", title: "remove from forcing.series", onclick: () => setOv("forcing", "series", asList(eff("forcing", "series", "")).filter(x => x !== name)) }, "×")), fields(g, ["method"])];
  if (derived) {
    const from = String(p.from || "");
    return h("div", { class: "fb-card" }, ...head,
      fields(g, ["from", ...(DERIVED_KEYS[from] || ["factor", "region", "anomaly"]), "units"]),
      note(from === "air_anomaly" ? `Set by the climate each update, after the atmosphere is composed, as ${fmtNum(num(p.factor))} × the mean over ${p.region === "all" ? "the climate grid" : `region "${p.region}"`} of the annual sea-level air temperature anomaly to the reference (no lag). No curve here; a term on it may touch to and so only.`
        : `The ${from === "obm_north" ? "northern" : "southern"} box temperature of the ocean box model${p.anomaly ? " minus its initial value" : ""}, set by the yelmox_bipolar bridge each step. No curve here; a term on it may touch to and so only.`));
  }
  const isfb = FEEDBACK.includes(method);
  const keys = isfb ? FEEDBACK_KEYS : (METHOD_KEYS[method] || ["f_min", "f_max"]);
  const opts = { series_file: { options: Object.keys(D.series).sort() } };
  if (method === "series" && p.series_file && D.series[p.series_file] && D.series[p.series_file].kind === "netcdf")
    opts.series_var = { options: Object.keys(D.series[p.series_file].vars) };
  const content = [...head, fields(g, keys, opts),
    details("series-extra-" + name, "noise, window, units, restarts", fields(g, ["sigma", "units", "restart_every_df", "with_kill", "win_t0", "win_t1", "win_taper", "win_fac", "win_ref", "win_above"]))];
  const s = num(p.df_sign) >= 0 ? 1 : -1;
  const v = [];
  if (isfb) v.push(`${method}: response-driven, no closed form. The curve is the fastest path, f_start + s·df_dt_max·τ_e from f_start = ${s > 0 ? "f_min" : "f_max"} after dt_init = ${fmtNum(p.dt_init)} yr; the controller slows it while the ice volume responds faster than tol${p.with_kill ? ", and the run stops once it equilibrates at a bound (with_kill)" : ""}.`);
  else if (method === "series") {
    v.push(`F: ${p.series_file}${p.series_var ? ":" + p.series_var : ""}`);
    if (p.norm) v.push(`normalised: t_0 = ${asList(p.norm_t0).map(fmtNum).join("…")}, t_1 = ${asList(p.norm_t1).map(fmtNum).join("…")} (a time, or the mean over a window)` + (ser && ser.norm ? ` → F̄[t_0] = ${fmtV(ser.norm[0])}, F̄[t_1] = ${fmtV(ser.norm[1])}` : ""));
    if (!ser) v.push("series not found in the exported data");
  } else if (METHOD_EQ[method]) {
    v.push(`τ_e = max(0, τ − τ(t_init) − dt_init); s = ${s}`);
    if ((method === "const" || method === "ramp-slope") && num(p.df_dt_max) > 0) v.push(`from ${fmtNum(s > 0 ? p.f_min : p.f_max)} to ${fmtNum(s > 0 ? p.f_max : p.f_min)} over ${fmtNum(Math.abs(num(p.f_max) - num(p.f_min)) / num(p.df_dt_max))} yr`);
    if (method === "sin") v.push("t_e: the time since the start of the run (held at 0 during dt_init)");
    if (tl.method === "const") v.push("the forcing time is constant: τ_e = 0, the series stays at f_start (tstep_method = rel advances it)");
  }
  if (num(p.win_t1) > num(p.win_t0)) v.push(`window ${fmtNum(p.win_t0)} < τ < ${fmtNum(p.win_t1)}: ${name} → ${fmtNum(p.win_ref)} + ${fmtNum(p.win_fac)}·(${name} − ${fmtNum(p.win_ref)})${p.win_above ? " where " + name + " > win_ref" : ""}`);
  if (num(p.sigma) > 0) v.push(`noise: + ε, ε ~ N(0, ${fmtNum(p.sigma)}²) every step`);
  if (v.length) content.push(h("ul", { class: "fb-list" }, ...v.map(x => h("li", {}, x))));
  if (f.some(Number.isFinite)) {
    const rules = (!isfb && method !== "series") ? [{ y: num(p.f_min), label: "f_min" }, { y: num(p.f_max), label: "f_max" }] : [];
    content.push(smallPlot([{ label: name, t: tl.t, f }], { y: name, rules }));
  }
  if (ser) {
    const vr = tl.method === "const" ? [{ x: tl.tc, label: "tstep_const" }] : [];
    if (p.norm) { for (const t of asList(p.norm_t0).map(num)) vr.push({ x: t, label: "t0", stroke: "#393" }); for (const t of asList(p.norm_t1).map(num)) vr.push({ x: t, label: "t1", stroke: "#393" }); }
    content.push(details("series-file-" + name, `the file (${ser.n} points; the run in yellow)`,
      smallPlot([{ label: "F", t: ser.time, f: ser.f }], { x: "forcing time τ [yr]", y: p.norm ? "normalised" : "F", vrules: vr, band: tl.method !== "const" ? { x0: tl.tfmin, x1: tl.tfmax } : null }),
      ser.header ? note(ser.header.slice(0, 160)) : null));
  }
  return h("div", { class: "fb-card" }, ...content);
}

function secSeries(tl) {
  const names = asList(eff("forcing", "series", ""));
  return section("series", "Scalar series", "named functions of τ: a file, a constant, a ramp, a controller, or derived",
    h("div", { class: "fb-card fb-card-flat" }, namesEditor("forcing", "series", named("series_"), "forcing.series"),
      details("var-defs", "where the named groups are read", note("[series_<name>], [term_<name>], [snapshot_<name>]: from the par file, then the files of forcing.var_defs (the forcing library). -p on a library group copies it into the run's par file."), fields("forcing", ["var_defs"]))),
    ...names.map(n => seriesCard(n, tl)));
}

// ------------------------------------------------- left: terms, snapshots

// One term [term_<name>] as a record.
function termOf(name) {
  const g = "term_" + name;
  return { name, g, series: String(eff(g, "series", "")), pattern: String(eff(g, "pattern", "snapshots")), snap0: String(eff(g, "snap0", "")), snap1: String(eff(g, "snap1", "")),
    fields: asList(eff(g, "fields", "")), gain: { tas: num(eff(g, "f_ta", 1)), to: num(eff(g, "f_to", 1)), so: num(eff(g, "f_so", 1)) } };
}
const boundFields = (snap) => FIELDS.filter(f => asList(eff("snapshot_" + snap, f, "")).length);
// time-dependent: a bound varslice group with a time axis (time = 1 t0 t1 dt ...)
const snapTdep = (snap) => FIELDS.some(f => asList(eff("snapshot_" + snap, f, "")).some(b => { const c = D.climdb[b]; return c && Array.isArray(c.time) && c.time[0] >= 1; }));
const snapSym = (f, snap) => `${FSYM[f]}_{${sub(snap)}}` + (snapTdep(snap) ? "(\\tau)" : "");
const wsym = (t) => t.series === "" ? "" : `${sym(t.series)}(\\tau)\\,`;
const gsym = (x) => x === 1 ? "" : (x === -1 ? "-" : `${fmtNum(x)}\\,`);
// The contribution of a term to a field, on its own (backends without a composition).
function termRhs(t, f) {
  if (t.pattern === "replace") return `${FSYM[f]} = ${sym(t.series)}(\\tau)\\quad\\text{at every depth, in place of the reference}`;
  return `${FSYM[f]} \\leftarrow ${FSYM[f]}` + termPart(t, f);
}
// The contribution of a term as a factor/summand of the composed equation.
function termPart(t, f) {
  if (t.pattern === "uniform") {
    if (f === "pr") return `\\,e^{\\beta_p\\,${t.gain.tas === 1 ? "" : "(" + fmtNum(t.gain.tas) + ")\\,"}${sym(t.series)}(\\tau)}`;
    const g = t.gain[f];
    return ` ${g < 0 ? "-" : "+"} ${gsym(Math.abs(g))}${sym(t.series)}(\\tau)`;
  }
  if (f === "pr") return t.snap0 ? `\\,\\Big[1 + ${wsym(t)}\\Big(\\frac{${snapSym(f, t.snap1)}}{${snapSym(f, t.snap0)}} - 1\\Big)\\Big]`
    : (t.series === "" ? `\\,${snapSym(f, t.snap1)}` : `\\,\\big[1 + ${wsym(t)}(${snapSym(f, t.snap1)} - 1)\\big]`);
  return t.snap0 ? ` + ${wsym(t)}\\big[${snapSym(f, t.snap1)} - ${snapSym(f, t.snap0)}\\big]` : ` + ${wsym(t)}${snapSym(f, t.snap1)}`;
}
const PATTERN_HELP = { snapshots: "the series × (snap1 − snap0)", uniform: "the series × a gain, the same everywhere", replace: "the field = the series, in place of the reference" };

function termWarnings(t, ctx) {
  const { snames, ref, withSnaps, cl } = ctx;
  const sinfo = t.series === "" ? null : (snames.includes(t.series) ? seriesOnRun(t.series, ctx.tl) : null);
  const derived = !!(sinfo && sinfo.derived);
  const w = [];
  if (!t.fields.length) w.push("no fields: the term does nothing");
  for (const f of t.fields) if (!FIELDS.includes(f)) w.push(`field "${f}" is not one of ${FIELDS.join(", ")}`);
  if (t.series !== "" && !snames.includes(t.series)) w.push(`series "${t.series}" is not in forcing.series`);
  if (derived && t.fields.some(f => f === "tas" || f === "pr")) w.push(`the derived series "${t.series}" is set after the atmosphere: this term may touch to and so only`);
  if (t.pattern === "replace") {
    if (t.series === "") w.push("replace needs a series");
    if (t.fields.some(f => f === "tas" || f === "pr")) w.push("replace: to and so only");
  }
  if (t.pattern === "snapshots") {
    if (!t.snap1) w.push("no snap1");
    for (const s of [t.snap0, t.snap1]) {
      if (!s) continue;
      if (!groupExists("snapshot_" + s)) { w.push(`snapshot "${s}": no group [snapshot_${s}] in the par file or the library`); continue; }
      for (const f of t.fields) if (FIELDS.includes(f) && !boundFields(s).includes(f)) w.push(`field ${f} is not bound by snapshot ${s}`);
    }
    if (derived) w.push("a derived series cannot drive a snapshot pattern (the snapshots are composed before it is set)");
  }
  if (withSnaps && groupExists("snapshot_" + ref)) for (const f of t.fields) if (FIELDS.includes(f) && !boundFields(ref).includes(f) && (cl === "snapshots" || f === "to" || f === "so")) w.push(`field ${f} is not bound by the reference snapshot ${ref}`);
  if (cl === "esm" && String(eff("comps", "surface_method")) === "climate" && t.fields.some(f => f === "tas" || f === "pr")) w.push("the ESM supplies the smb (surface_method = climate): terms on tas/pr are refused");
  if (cl === "rembo" && t.fields.some(f => f === "tas" || f === "pr") && t.pattern !== "uniform") w.push("rembo takes uniform tas terms only (its summer anomaly); snapshot terms on the atmosphere are not applied");
  return { w, derived };
}

function termCard(name, ctx) {
  const g = "term_" + name, src = groupSource(g);
  if (!src) return h("div", { class: "fb-card" }, warn(`[${g}]: no such group and no term template`));
  const t = termOf(name);
  const { w, derived } = termWarnings(t, ctx);
  const content = [h("div", { class: "fb-cardhead" }, h("code", {}, `[${g}]`), h("span", { class: "fb-src fb-src-" + src }, src === "new" ? "new (template term)" : src),
    h("span", { class: "fb-kind" }, t.pattern), h("span", { class: "fb-help" }, PATTERN_HELP[t.pattern] || ""),
    h("button", { class: "fb-x", title: "remove from forcing.terms", onclick: () => setOv("forcing", "terms", asList(eff("forcing", "terms", "")).filter(x => x !== name)) }, "×"))];
  content.push(fields(g, ["series", "pattern"], { series: { options: ["", ...ctx.snames], labels: { "": "\"\" (the constant 1)" } } }));
  const marks = {};
  if (t.pattern === "snapshots") {
    content.push(fields(g, ["snap0", "snap1"], { snap0: { options: ["", ...ctx.candSnap], labels: { "": "\"\" (none: snap1 is the anomaly or ratio)" } }, snap1: { options: ctx.candSnap } }));
    for (const f of FIELDS) {
      const unbound = [t.snap0, t.snap1].filter(s => s && groupExists("snapshot_" + s) && !boundFields(s).includes(f));
      if (unbound.length) marks[f] = `not bound by ${unbound.join(", ")}`;
    }
  }
  if (t.pattern === "replace") for (const f of ["tas", "pr"]) marks[f] = "replace: to and so only";
  if (derived) for (const f of ["tas", "pr"]) marks[f] = `derived series ${t.series}: to and so only`;
  content.push(checkList(g, "fields", FIELDS, marks));
  if (t.pattern === "uniform") {
    const gk = ["tas", "to", "so"].filter(f => t.fields.includes(f)).map(f => GAIN[f]);
    if (gk.length) content.push(fields(g, gk));
    if (t.fields.includes("pr")) content.push(note(`pr follows the tas gain: × exp(β_p f_ta ${t.series}), β_p = snapshots.f_p (fixed rule); leave pr out for no precipitation response.`));
  }
  content.push(...w.map(warn));
  return h("div", { class: "fb-card" }, ...content);
}

function snapshotCard(name, role) {
  const g = "snapshot_" + name, src = groupSource(g);
  if (!groupExists(g)) return h("div", { class: "fb-card" }, warn(`[${g}] (${role}): no such group in the par file or the library`));
  const bound = ["tas", "pr", "zs", "to", "so"].map(fld => {
    const v = asList(eff(g, fld, ""));
    return v.length ? h("li", {}, h("code", {}, fld), ": ", v.map(b => { const c = D.climdb[b]; return c ? `${b} (${c.filename.replace(/^fesmdata\//, "")}${Array.isArray(c.time) && c.time[0] >= 1 ? `, in time ${fmtNum(c.time[1])} to ${fmtNum(c.time[2])}` : ""})` : b + " (not in input/clim_db.nml)"; }).join("; ")) : null;
  });
  return card(`[${g}]`, src, h("span", { class: "fb-help" }, role + (snapTdep(name) ? "; depends on τ" : "")),
    details("snap-" + name, "fields bound", h("ul", { class: "fb-list" }, ...bound), fields(g, ["exp", "monthly", "clim_sigma", "ocn_sigma"])));
}

function termCtx(tl) {
  const cl = climate();
  return { tl, cl, withSnaps: ["snapshots", "rembo"].includes(cl), snames: asList(eff("forcing", "series", "")), tnames: asList(eff("forcing", "terms", "")),
    ref: String(eff("forcing", "reference", "pd")), candSnap: named("snapshot_") };
}
const runTerms = (ctx) => ctx.tnames.filter(n => groupSource("term_" + n)).map(termOf);
// The snapshots of the run: the reference and the snap0/snap1 of the terms.
function runSnaps(ctx) {
  const snaps = [ctx.ref];
  for (const t of runTerms(ctx)) if (t.pattern === "snapshots") for (const s of [t.snap0, t.snap1]) if (s && !snaps.includes(s)) snaps.push(s);
  return snaps;
}

function secTerms(ctx) {
  const { cl, withSnaps, ref, candSnap } = ctx;
  return section("terms", "Terms", "series × pattern on explicit fields, added to the reference",
    h("div", { class: "fb-card fb-card-flat" },
      note(cl === "snapshots" ? "tas, to, so add, pr multiplies (fixed by field). Applied in this order: the snapshot terms, the uniform terms on tas/pr, the derived series, the ocean terms."
        : cl === "rembo" ? "REMBO supplies the atmosphere and takes the uniform tas terms as its summer anomaly (ΔT_sum = Σ f_ta α, ΔT_ann = 1.3 ΔT_sum); the snapshot backend supplies the ocean from the reference snapshot plus the ocean terms."
        : "The ESM composes the atmosphere and the ocean from its own anomalies; the terms add to its result: uniform terms on tas/pr unless the ESM supplies the smb, ocean terms at the shelf base."),
      withSnaps ? fields("forcing", ["reference"], { reference: { options: candSnap } }) : null,
      namesEditor("forcing", "terms", named("term_"), "forcing.terms")),
    ...ctx.tnames.map(n => termCard(n, ctx)));
}

function secSnapshots(ctx) {
  if (!ctx.withSnaps) return null;
  const terms = runTerms(ctx);
  const roles = (s) => [s === ctx.ref ? "the reference" : null, ...terms.filter(t => t.pattern === "snapshots" && (t.snap0 === s || t.snap1 === s)).map(t => `${t.snap0 === s ? "snap0" : "snap1"} of ${t.name}`)].filter(Boolean).join(", ");
  return section("snapshots", "Snapshots", "the reference and the snap0/snap1 of the terms",
    h("div", { class: "fb-card fb-card-flat" }, fields("snapshots", ["lapse", "f_p"]),
      details("snapshots-extra", "combine, grid_global, fields, f_p_ne", fields("snapshots", ["combine", "grid_global", "fields", "f_p_ne"]))),
    ...runSnaps(ctx).map(s => snapshotCard(s, roles(s))));
}

function secEsm(tl) {
  if (climate() !== "esm") return null;
  const T = tl.T;
  const periods = ["time_ref", "time_esm_ref", "time_hist", "time_proj"].map(k => { const v = asList(eff(T, k, "")).map(num); return v.length === 2 ? { k, x0: v[0], x1: v[1] } : null; }).filter(Boolean);
  if (tl.cal) periods.push({ k: "run (" + T + ")", x0: tl.method === "const" ? tl.tc : tl.t0, x1: tl.method === "const" ? tl.tc : tl.t1 });
  const fig = Plot.plot({ width: WL, height: 150, marginLeft: 110, x: { label: "year CE", grid: true }, y: { label: null },
    marks: [Plot.barX(periods, { x1: "x0", x2: "x1", y: "k", fill: "#69b" }), Plot.ruleX([tl.method === "const" ? tl.tc : NaN], { stroke: "#c33", strokeDasharray: "4 3" })] });
  return section("esm", "ESM forcing", "a reference climatology plus the anomalies of an ESM",
    h("div", { class: "fb-card" },
      fields("esm", ["experiment", "esm_name", "use_esm", "use_hist", "use_proj", "use_var", "basins", "par_file"]),
      fields("esm", ["lapse", "f_p", "f_ocn", "f_polar", "dT_threshold"]),
      h("figure", { class: "fb-fig" }, fig)));
}

// ----------------------------------------- left: transforms, sea level, orbit

function secTransforms() {
  return section("transforms", "Transforms", "after the composition; no series enter here",
    card("[taper]", null, h("span", { class: "fb-help" }, "the atmosphere towards the reference with a weight field"), fields("taper", ["active", "weight", "outside"])),
    card("[bias]", null, h("span", { class: "fb-help" }, "an additive temperature field and a precipitation factor field"), fields("bias", ["active", "dtas", "fpr", "f_tas", "f_pr"])));
}

function secExtras() {
  const snames = asList(eff("forcing", "series", ""));
  const content = [];
  const m = String(eff("barysealevel", "method", "fastiso"));
  const sl = [fields("barysealevel", ["method", "bsl_init", "sl_path", "sl_name"], { sl_path: { options: ["none", ...Object.keys(D.series).sort()] } })];
  if (m === "file") {
    const path = String(eff("barysealevel", "sl_path", "")), vname = String(eff("barysealevel", "sl_name", ""));
    if (!seriesOf(path, vname)) sl.push(warn(`${path}${vname && vname !== "none" ? ":" + vname : ""}: not in the exported series (netCDF: set sl_name)`));
  }
  content.push(card("[barysealevel]", null, h("span", { class: "fb-help" }, "sea level"), ...sl));
  const om = String(eff("orbit", "method", "fixed_time"));
  content.push(card("[orbit]", null, h("span", { class: "fb-help" }, "insolation of the surface models"), fields("orbit", ["method", ...(om === "fixed_time" ? ["time"] : om === "fixed_orbit" ? ["ecc", "obl", "per"] : [])])));
  if (groupSource("ice_bnd")) {
    const im = String(eff("ice_bnd", "method"));
    const ib = [fields("ice_bnd", ["active", "method", "H_ice", ...(im === "scaled" ? ["index", "slice_time"] : [])], { index: { options: ["", ...snames] } })];
    if (eff("ice_bnd", "active", false) && im === "scaled" && !snames.includes(String(eff("ice_bnd", "index", "")))) ib.push(warn(`ice_bnd.index = "${eff("ice_bnd", "index", "")}" is not in forcing.series`));
    content.push(card("[ice_bnd]", null, h("span", { class: "fb-help" }, "prescribed ice in the fixed and relaxation cells"), ...ib));
  }
  return section("extras", "Also at the forcing time", "sea level, orbit, prescribed ice", ...content);
}

// ------------------------------------------------- right: the equation chain

function rightEquations(tl, ctx) {
  const t = evalTime(tl), tau = tauAt(tl, t);
  const out = [];
  const { cl, withSnaps, snames, ref } = ctx;
  const terms = runTerms(ctx);
  const sv = {};   // series values at t
  for (const n of snames) { const s = seriesOnRun(n, tl); sv[n] = s.derived ? null : atT(tl, s.f, t); }
  const tauS = `${fmtNum(+tau.toPrecision(6))} yr`;

  // 0. forcing time, with the evaluation slider
  const timeEq = tl.method === "const" ? `\\tau = t_{\\mathrm{const}} = ${fmtK(tl.tc)}` : tl.method === "rel" ? `\\tau = t` : `\\tau = t_{\\mathrm{cal}} - ${fmtNum(tl.cref)}`;
  const lo = Math.min(tl.t0, tl.t1), hi = Math.max(tl.t0, tl.t1);
  const slider = h("input", { type: "range", min: lo, max: hi, step: (hi - lo) / 400 || 1, value: t, oninput: (e) => { S.tsel = Number(e.target.value); renderRight(); } });
  out.push(eqblock("0 · Forcing time",
    tex(`${timeEq},\\qquad t \\in [${fmtK(tl.t0)},\\ ${fmtK(tl.t1)}]\\ \\mathrm{yr}\\quad(\\texttt{${tl.method}})`),
    h("div", { class: "fb-taurow" }, slider, h("span", { class: "fb-tauval" }, `t = ${fmtNum(+t.toPrecision(6))} yr → τ = ${tauS}`)),
    vals("Every series, time-dependent snapshot, sea-level file and the orbit (forcing_time) read τ. The values below are at this τ.")));

  // 1. series
  if (snames.length) {
    const b = eqblock("1 · Scalar series");
    for (const n of snames) {
      if (!groupSource("series_" + n)) { b.append(warn(`series "${n}": no group [series_${n}]`)); continue; }
      const p = seriesPar("series_" + n);
      b.append(tex(seriesEq(n, p) + (sv[n] === null ? "" : `\\qquad ${sym(n)}(${fmtK(+tau.toPrecision(6))}) = ${fmtV(sv[n])}`)));
    }
    out.push(b);
  }

  // 2. snapshots, 3. composition, 4. to the surface
  if (withSnaps) {
    const lapse = asList(eff("snapshots", "lapse", [0.008, 0.0065]));
    out.push(eqblock("2 · Snapshots, reduced to sea level once",
      tex(`T^{\\mathrm{sl}}_j = T_j + \\Gamma_m z_j,\\qquad P^{\\mathrm{sl}}_j = P_j\\,e^{\\beta_p \\Gamma_m z_j},\\qquad \\Theta_j(z),\\ S_j(z)\\ \\text{as given}`),
      vals(`Γ = ${fmtNum(lapse[0])}, ${fmtNum(lapse[1])} K/m (annual, summer; monthly cosine between them), β_p = ${fmtNum(num(eff("snapshots", "f_p", 0.05)))} /K; j ∈ {${runSnaps(ctx).map(s => s + (snapTdep(s) ? "(τ)" : "")).join(", ")}}, reference ${ref}`)));
    const b3 = eqblock(cl === "rembo" ? "3 · Composition (REMBO: the atmosphere; the snapshots: the ocean)" : "3 · Composition at sea level: reference + Σ series × pattern");
    if (!groupExists("snapshot_" + ref)) b3.append(warn(`the reference snapshot "${ref}" has no group [snapshot_${ref}]`));
    else {
      const refF = boundFields(ref);
      for (const f of FIELDS) {
        if (!refF.includes(f)) continue;
        if (cl === "rembo" && (f === "tas" || f === "pr")) continue;
        const tf = terms.filter(t => t.fields.includes(f));
        const rp = tf.find(t => t.pattern === "replace");
        let e = `${FSYM[f]}(${f === "tas" || f === "pr" ? "x,m" : "x,z"},\\tau) = ` + (rp ? `${sym(rp.series)}(\\tau)` : snapSym(f, ref));
        for (const t of tf) if (t.pattern !== "replace") e += termPart(t, f);
        b3.append(tex(e));
      }
      if (cl === "rembo") { const u = terms.filter(t => t.pattern === "uniform" && t.fields.includes("tas")); if (u.length) b3.append(tex(`\\Delta T_{\\mathrm{sum}} = ${u.map(t => gsym(t.gain.tas) + sym(t.series) + "(\\tau)").join(" + ")},\\qquad \\Delta T_{\\mathrm{ann}} = 1.3\\,\\Delta T_{\\mathrm{sum}}`)); }
      for (const f of ["to", "so"]) if (terms.filter(t => t.pattern === "replace" && t.fields.includes(f)).length > 1) b3.append(warn(`more than one replace term on ${f}`));
    }
    const num3 = terms.map(t => t.series === "" ? `${t.name}: weight 1` : sv[t.series] === null ? `${t.name}: ${t.series} set at run time` : `${t.name}: ${t.series} = ${fmtV(sv[t.series])}` + (t.pattern === "uniform" && t.fields.includes("tas") ? ` → ΔT = ${fmtV(t.gain.tas * sv[t.series])} K` : ""));
    if (num3.length) b3.append(vals(`at τ = ${tauS}:  ` + num3.join("  ·  ")));
    out.push(b3);
    if (cl === "snapshots") out.push(eqblock("4 · To the model surface z_s(x,t)", tex(`T(x,m) = T^{\\mathrm{sl}} - \\Gamma_m z_s,\\qquad P(x,m) = P^{\\mathrm{sl}}\\,e^{-\\beta_p \\Gamma_m z_s}`)));
  } else if (cl === "esm") {
    const T = tl.T, useEsm = eff("esm", "use_esm", false), exp = String(eff("esm", "experiment", "ctrl"));
    const b = eqblock("2 · The ESM climate");
    if (useEsm && exp !== "ctrl") b.append(tex(`X(t) = X_{\\mathrm{ref}} + \\big(X_{\\mathrm{ESM}}(t) - \\overline{X_{\\mathrm{ESM}}}^{\\,\\mathrm{time\\_esm\\_ref}}\\big),\\quad \\text{historical to } ${fmtNum(asList(eff(T, "time_hist", [0, 0]))[1])},\\ \\text{projection from } ${fmtNum(asList(eff(T, "time_proj", [0, 0]))[0])}`),
      vals("Air temperature and smb anomalies are differences, precipitation a ratio, the ocean annual means at the shelf base."));
    else if (!useEsm && exp === "1pctCO2") b.append(tex(`\\Delta T_a(t) = f_{\\mathrm{polar}}\\,\\min\\!\\big(17\\,\\mathrm{K}\\,\\tfrac{t - 1995}{305},\\ \\Delta T_{\\mathrm{thr}}\\big),\\qquad \\Delta T_o = f_{\\mathrm{ocn}}\\,\\Delta T_a`),
      vals(`f_polar = ${fmtNum(eff("esm", "f_polar"))}, ΔT_thr = ${fmtNum(eff("esm", "dT_threshold"))}, f_ocn = ${fmtNum(eff("esm", "f_ocn"))}: homogeneous warming (use_esm = False, experiment 1pctCO2)`));
    else b.append(vals(`experiment = "${exp}", use_esm = ${useEsm}: the reference climatology only (spin-up).`));
    out.push(b);
    if (terms.length) { const b3 = eqblock("3 · Terms on the ESM result"); for (const t of terms) for (const f of t.fields) if (FIELDS.includes(f)) b3.append(tex(`\\text{${t.name}:}\\quad ` + termRhs(t, f))); out.push(b3); }
  }

  // 5. transforms
  const tr = [];
  if (eff("taper", "active", false)) tr.push(tex(`X \\leftarrow X_0(z_s) + w(x)\\,\\big(X - X_0(z_s)\\big)` + (String(eff("taper", "outside")) === "mean" ? ` + (1-w)\\,\\overline{(X - X_0)}^{\\,w=1}` : "") + `,\\qquad X \\in \\{T, P\\}`),
    vals(`taper: w = ${eff("taper", "weight", "")} (1: the composed climate, 0: the reference at the current surface)`));
  if (eff("bias", "active", false)) {
    const tdep = climTdep(eff("bias", "dtas", "")) || climTdep(eff("bias", "fpr", ""));
    const a = tdep ? "(x,m,\\tau)" : "(x,m)";
    tr.push(tex(`T \\leftarrow T + ${gsym(num(eff("bias", "f_tas")))}b_T${a},\\qquad P \\leftarrow P\\,\\big(1 + ${gsym(num(eff("bias", "f_pr")))}(b_P${a} - 1)\\big)`), vals(`bias: b_T = ${eff("bias", "dtas", "")}, b_P = ${eff("bias", "fpr", "")}` + (tdep ? " (a time axis: sliced at τ)" : "")));
  }
  out.push(eqblock("5 · Transforms", ...(tr.length ? tr : [vals("none active (taper, bias)")])));

  // 6. also at the forcing time
  const b6 = eqblock("6 · Also at the forcing time");
  const om = String(eff("orbit", "method", "fixed_time"));
  b6.append(om === "forcing_time" ? tex(`I = I\\big(\\varphi, d;\\ \\mathrm{orbit}(\\tau)\\big)\\quad\\text{Laskar (2004)}`)
    : om === "fixed_time" ? tex(`I = I\\big(\\varphi, d;\\ \\mathrm{orbit}(${fmtK(num(eff("orbit", "time")))}\\ \\mathrm{yr})\\big)\\quad\\text{fixed, Laskar (2004)}`)
    : om === "fixed_orbit" ? tex(`I = I(\\varphi, d;\\ e, \\varepsilon, \\varpi),\\qquad e = ${fmtNum(eff("orbit", "ecc"))},\\ \\varepsilon = ${fmtNum(eff("orbit", "obl"))}^\\circ,\\ \\varpi = ${fmtNum(eff("orbit", "per"))}^\\circ`)
    : warn(`orbit.method = "${om}": forcing_time, fixed_time or fixed_orbit`));
  const m = String(eff("barysealevel", "method", "fastiso"));
  if (m === "const") b6.append(tex(`\\eta = ${fmtNum(eff("barysealevel", "bsl_init"))}\\ \\mathrm{m}`));
  else if (m === "fastiso") b6.append(tex(`\\eta(t) = \\eta_0 - \\Delta V_{\\mathrm{ice}}(t)/A_{\\mathrm{ocean}}\\quad\\text{(the FastIsostasy domains)}`));
  else if (m === "file") { const ser = seriesOf(String(eff("barysealevel", "sl_path", "")), String(eff("barysealevel", "sl_name", ""))); b6.append(tex(`\\eta = \\eta_{\\mathrm{file}}(\\tau)` + (ser ? ` = ${fmtV(interp(ser, tau))}\\ \\mathrm{m}` : ""))); }
  if (groupSource("ice_bnd") && eff("ice_bnd", "active", false)) {
    const im = String(eff("ice_bnd", "method")), s = String(eff("ice_bnd", "index", ""));
    b6.append(im === "scaled" ? tex(`H_b(x) = ${s ? sym(s) + "(\\tau)" : "?"}\\,H(x, ${fmtK(num(eff("ice_bnd", "slice_time")))})`) : tex(`H_b(x) = H(x, \\tau)`));
  }
  out.push(b6);
  return out;
}

// ------------------------------------------------------ right: over the run

function rightRun(tl, ctx) {
  const t = evalTime(tl);
  const terms = runTerms(ctx);
  const content = [];
  const ser = ctx.snames.map(n => ({ n, s: seriesOnRun(n, tl) })).filter(x => x.s.f && x.s.f.some(Number.isFinite)).map(x => ({ label: x.n, t: tl.t, f: x.s.f }));
  if (ser.length) content.push(linePlot(ser, { y: "series", tmark: t }));
  const uni = terms.filter(x => x.pattern === "uniform");
  if (uni.length) {
    const sums = {};
    for (const x of uni) {
      const w = termWeight(x, tl, ctx.snames);
      if (!w) continue;
      for (const [f, lab] of [["tas", ctx.cl === "rembo" ? "ΔT_sum [K]" : "ΔT_a [K]"], ["to", "ΔΘ [K]"], ["so", "ΔS [psu]"]]) if (x.fields.includes(f)) sums[lab] = (sums[lab] || tl.t.map(() => 0)).map((y, i) => y + x.gain[f] * w[i]);
    }
    const s = Object.entries(sums).map(([label, f]) => ({ label, t: tl.t, f }));
    if (s.length) content.push(linePlot(s, { y: "uniform anomalies", height: 160, tmark: t }));
  }
  if (String(eff("barysealevel", "method")) === "file") {
    const sl = seriesOf(String(eff("barysealevel", "sl_path", "")), String(eff("barysealevel", "sl_name", "")));
    if (sl) content.push(linePlot([{ label: "η [m]", t: tl.t, f: tl.tf.map(x => interp(sl, x)) }], { y: "sea level [m]", height: 140, tmark: t }));
  }
  if (!content.length) content.push(vals(tl.method === "const" ? "The forcing time is constant: nothing changes over the run." : "No tabulated series in the run (derived series are set by the model)."));
  return eqblock("Over the run (model time; the line: the evaluation time)", ...content);
}

// --------------------------------------------- right: the run as text (tabs)

// The named groups the configuration references: the series, the terms, the
// reference and the snapshots of the terms (in this order, after the base's groups).
function referencedGroups() {
  const out = [];
  const add = (g) => { if (!out.includes(g)) out.push(g); };
  for (const n of asList(eff("forcing", "series", ""))) add("series_" + n);
  for (const n of asList(eff("forcing", "terms", ""))) add("term_" + n);
  if (["snapshots", "rembo"].includes(climate())) {
    add("snapshot_" + String(eff("forcing", "reference", "pd")));
    for (const n of asList(eff("forcing", "terms", ""))) {
      if (!groupSource("term_" + n)) continue;
      const t = termOf(n);
      if (t.pattern === "snapshots") for (const s of [t.snap0, t.snap1]) if (s) add("snapshot_" + s);
    }
  }
  return out;
}
function fmtP(v) {
  if (typeof v === "boolean") return v ? "True" : "False";
  if (typeof v === "number") return fmtNum(v);
  const s = String(v);
  return /^[\w.\/{}:~&|,-]*$/.test(s) ? s : JSON.stringify(s);
}
function fmtNml(v) {
  if (Array.isArray(v)) return v.length ? v.map(fmtNml).join(" ") : "\"\"";
  if (typeof v === "boolean") return v ? "True" : "False";
  if (typeof v === "number") return fmtNum(v);
  return JSON.stringify(String(v));
}
// The groups of the run in par-file order: the base's, then the referenced ones.
function runGroups() {
  const order = [...base().order];
  const add = (g) => { if (!order.includes(g)) order.push(g); };
  referencedGroups().forEach(add);
  for (const id of S.ov.keys()) add(id.split(".")[0]);
  return order.filter(g => groupSource(g));
}
function overrides() {
  const out = [], blocks = [], problems = [];
  for (const g of runGroups()) {
    const ids = [...S.ov.keys()].filter(id => id.split(".")[0] === g);
    if (!ids.length) continue;
    const src = groupSource(g);
    if (src === "new" && !D.defaults[g]) {
      // runme -p cannot create a templated group: a namelist block for the par file
      const lines = ids.map(id => { const k = id.split(".")[1]; return `    ${k.padEnd(16)} = ${fmtNml(S.ov.get(id))}`; });
      blocks.push(`&${g}\n${lines.join("\n")}\n/`);
      problems.push(`[${g}] is new: -p cannot create it. Add the block below to the par file or to input/forcing_db.nml.`);
      continue;
    }
    for (const id of ids) {
      const k = id.split(".")[1], v = S.ov.get(id);
      out.push(Array.isArray(v) ? `${g}.${k}==${v.map(fmtP).join(",")}` : `${g}.${k}=${fmtP(v)}`);
    }
    if (src === "library") problems.push(`[${g}] comes from the library: -p on it needs runme ≥ 0.5.26 (it copies the group into the staged par file).`);
  }
  return { out: [...out, ...S.extra], blocks, problems };
}
function parseP(text) {
  // runme -p tokens (whitespace separated; quoted values allowed) -> overrides;
  // keys outside the exported groups (yelmo physics, regions, ...) are kept as
  // they are and passed through to the command
  const toks = text.match(/(?:[^\s"']+|"[^"]*"|'[^']*')+/g) || [];
  const errs = [];
  for (let tok of toks) {
    if (tok === "-p" || tok === "\\") continue;
    const m = tok.match(/^([\w]+)\.([\w]+)(==|=)(.*)$/s);
    if (!m) { errs.push(`not a key=value: ${tok}`); continue; }
    const [, g, k, op, raw] = m;
    const info = keyInfo(g, k);
    if (!info) { if (!S.extra.includes(tok)) S.extra.push(tok); continue; }
    const items = (op === "==" ? raw.split(",") : [raw]).map(s => s.replace(/^["']|["']$/g, ""));
    const refel = Array.isArray(info.ref) ? (info.ref[0] ?? "") : info.ref;
    const conv = (s) => typeof refel === "boolean" ? /^(true|t|1|\.true\.)$/i.test(s) : typeof refel === "number" ? num(s) : s;
    const v = (op === "==" || Array.isArray(info.ref)) ? items.map(conv) : conv(items[0]);
    if (sameValue(v, info.ref)) S.ov.delete(g + "." + k); else S.ov.set(g + "." + k, v);
  }
  return errs;
}

// Tab "par file": the forcing and run-control groups of the run, every key
// resolved, with where its value comes from.
function parFile() {
  const lines = [];
  const esc = (s) => s.replace(/&/g, "&amp;").replace(/</g, "&lt;");
  for (const g of runGroups()) {
    const src = groupSource(g);
    lines.push(`<span class="g">&amp;${g}</span>` + (src === "library" ? `   <span class="c">! from the forcing library (${esc(D.library.files.join(", "))})</span>` : src === "new" ? `   <span class="c">! new, from the template</span>` : ""));
    const keys = schemaKeys(g);
    const w = Math.max(...keys.map(k => k.length), 4);
    for (const k of keys) {
      const i = keyInfo(g, k);
      const c = i.src === "override" ? `! -p (${i.refsrc}: ${fmtNml(i.ref)})` : i.src === "default" ? "! default" : "";
      lines.push(`<span${i.src === "override" ? ' class="ov"' : ""}>    ${k.padEnd(w)} = ${esc(fmtNml(i.v))}${c ? `   <span class="c">${esc(c)}</span>` : ""}</span>`);
    }
    lines.push("/", "");
  }
  if (S.extra.length) lines.push(`<span class="c">! -p outside the forcing and run-control groups, passed through:</span>`, ...S.extra.map(x => `<span class="c">!   ${esc(x)}</span>`));
  return lines.join("\n");
}

// Tab "equation as text": the composition written as one line per field (a
// design idea of docs/_forcing-terms.md, not read by the model).
function eqText(ctx) {
  const terms = runTerms(ctx);
  const w = (t) => t.series === "" ? "" : t.series + "*";
  const rows = [];
  if (ctx.withSnaps && groupExists("snapshot_" + ctx.ref)) {
    const refF = boundFields(ctx.ref);
    for (const f of FIELDS) {
      if (!refF.includes(f)) continue;
      const ts = terms.filter(t => t.fields.includes(f));
      const rp = ts.find(t => t.pattern === "replace");
      let e = rp ? rp.series : ctx.ref;
      for (const t of ts) {
        if (t.pattern === "replace") continue;
        if (f === "pr") e += t.pattern === "snapshots" ? (t.snap0 ? ` * (1 + ${w(t)}(${t.snap1}/${t.snap0} - 1))` : ` * (1 + ${w(t)}(${t.snap1} - 1))`) : ` * exp(beta_p*${t.gain.tas === 1 ? "" : fmtNum(t.gain.tas) + "*"}${t.series})`;
        else if (t.pattern === "snapshots") e += t.snap0 ? ` + ${w(t)}(${t.snap1} - ${t.snap0})` : ` + ${w(t)}${t.snap1}`;
        else { const g = t.gain[f]; e += ` ${g < 0 ? "-" : "+"} ${Math.abs(g) === 1 ? "" : fmtNum(Math.abs(g)) + "*"}${t.series}`; }
      }
      rows.push([f, e]);
    }
  } else for (const t of terms) for (const f of t.fields) rows.push([f, `${f} ${t.pattern === "replace" ? "= " + t.series : (t.pattern === "uniform" ? "+ " + (t.gain[f] === 1 ? "" : fmtNum(t.gain[f]) + "*") + t.series : "+ " + w(t) + "(" + t.snap1 + (t.snap0 ? " - " + t.snap0 : "") + ")")}`]);
  const l = ["&forcing", `    series   = ${fmtNml(ctx.snames)}`];
  for (const [f, e] of rows) l.push(`    ${f.padEnd(8)} = "${e}"`);
  l.push(`    var_defs = ${fmtNml(asList(eff("forcing", "var_defs", "")))}`, "/", "");
  for (const n of ctx.snames) {
    const g = "series_" + n;
    if (!groupSource(g)) continue;
    l.push(`&${g}`);
    for (const k of schemaKeys(g)) { const i = keyInfo(g, k); if (i.src !== "default") l.push(`    ${k} = ${fmtNml(i.v)}`); }
    l.push("/", "");
  }
  return l.join("\n").replace(/&/g, "&amp;").replace(/</g, "&lt;");
}

function rightOutput(ctx) {
  const tabs = [["par", "par file"], ["eq", "equation as text"], ["p", "-p and the runme command"]];
  const head = h("div", { class: "fb-tabs" }, ...tabs.map(([k, lab]) => h("button", { class: S.tab === k ? "on" : "", onclick: () => { S.tab = k; renderRight(); } }, lab)));
  const body = [];
  if (S.tab === "par") {
    body.push(note("The forcing and run-control groups of the run as the model reads them: the base par file, the library groups it names, every key resolved. Highlighted: the -p overrides."),
      h("pre", { class: "fb-cfg", html: parFile() }));
  } else if (S.tab === "eq") {
    body.push(note("A design idea, not implemented: the same composition written as one equation per field, the series by name (docs/_forcing-terms.md). The model reads the groups of the par-file tab. Open point: the pr line repeats what the tas line implies (ratio form, e^{β_p ΔT})."),
      h("pre", { class: "fb-cfg", html: eqText(ctx) }));
  } else {
    const { out, blocks, problems } = overrides();
    const outdir = S.outdir || `output/${S.preset || S.base.replace(/^yelmox_/, "").toLowerCase() + "-run"}`;
    const cmd = [`runme ${S.runopts} -e yelmox`, `-n ${base().file}`, `-o ${outdir}`, ...(out.length ? ["-p " + out.join(" \\\n      ")] : [])].join(" \\\n   ");
    const ta = h("textarea", { class: "fb-ta", rows: 3, placeholder: "paste a -p list, e.g. ctrl.tstep_method=rel forcing.series==orbital,dTo,ramp forcing.terms==orbital,ocean,ramp" });
    body.push(note("Only keys that differ from the base (or the library group, or the schema default); a library group is used by name. Keys outside the forcing groups (from a preset or a pasted list) are passed through."),
      h("div", { class: "fb-fields" },
        h("div", { class: "fb-field" }, h("label", {}, h("code", {}, "runme options")), h("input", { type: "text", class: "fb-text", value: S.runopts, onchange: (e) => { S.runopts = e.target.value; renderRight(); } }), h("span"), h("span", { class: "fb-help" }, "-r runs, -s stages, -rs submits; queue options (-q, -w, --omp)")),
        h("div", { class: "fb-field" }, h("label", {}, h("code", {}, "-o")), h("input", { type: "text", class: "fb-text", value: outdir, onchange: (e) => { S.outdir = e.target.value; renderRight(); } }), h("span"), h("span", { class: "fb-help" }, "the run folder"))),
      h("pre", { class: "fb-cfg" }, h("code", {}, cmd)),
      h("div", {}, h("button", { class: "fb-btn", onclick: () => navigator.clipboard.writeText(cmd.replace(/ \\\n\s+/g, " ")) }, "copy (one line)")),
      ...problems.map(warn),
      ...(blocks.length ? [h("pre", { class: "fb-cfg" }, h("code", {}, blocks.join("\n\n")))] : []),
      h("h4", {}, "Load a -p list"), ta,
      h("div", {}, h("button", { class: "fb-btn", onclick: () => { const errs = parseP(ta.value); S.preset = null; render(); if (errs.length) alert(errs.join("\n")); } }, "load -p list")));
  }
  return h("div", { class: "fb-eqblock fb-output" }, h("h3", {}, "The run, written out"), head, ...body);
}

const climTdep = (g) => { const c = D.climdb[String(g)]; return !!(c && Array.isArray(c.time) && c.time[0] >= 1); };

// --------------------------------------------------------------------- render

let LEFT = null, RIGHT = null;
function renderRight() {
  const tl = timeline(), ctx = termCtx(tl);
  const y = RIGHT.scrollTop;
  RIGHT.replaceChildren(...rightEquations(tl, ctx), rightRun(tl, ctx), rightOutput(ctx));
  RIGHT.scrollTop = y;
}
function render() {
  const y = window.scrollY;
  const tl = timeline(), ctx = termCtx(tl);
  root.replaceChildren();
  LEFT = h("div", { class: "fb-left" }, ...[secRun(tl), secSeries(tl), secTerms(ctx), secSnapshots(ctx), secEsm(tl), secTransforms(), secExtras()].filter(Boolean),
    h("p", { class: "fb-help" }, `Data: forcing.json generated ${D.generated}${D.commit ? " from " + D.commit : ""} (scripts/forcing-builder/export.jl).`));
  RIGHT = RIGHT || h("div", { class: "fb-right" });
  root.append(topBar(), h("div", { class: "fb-cols" }, LEFT, RIGHT));
  renderRight();
  window.scrollTo(0, y);
}

fetch(new URL("forcing-builder/forcing.json", document.baseURI)).then(r => r.json()).then(d => {
  D = d;
  TEMPLATES = (D.templates || []).map(t => [new RegExp(t.pattern), t.template]);
  const q = new URLSearchParams(location.search);
  const pre = (D.presets || []).find(p => p.name === q.get("preset"));
  if (pre) { loadPreset(pre); return; }
  S.base = D.bases[q.get("base")] ? q.get("base") : (D.bases.yelmox_Greenland ? "yelmox_Greenland" : Object.keys(D.bases)[0]);
  if (q.get("p")) parseP(q.get("p"));
  render();
}).catch(err => { root.replaceChildren(h("p", { class: "fb-warn" }, "forcing.json could not be loaded: " + err)); });
