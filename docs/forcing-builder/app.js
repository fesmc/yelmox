// Forcing builder: a view on one run's configuration (a base par file plus
// runme -p overrides, resolved against the forcing library and the schema
// defaults), with the forcing it imposes drawn on the run's timeline, and the
// runme command that reproduces it. Data: forcing.json (scripts/forcing-builder/
// export.jl). The index formulas mirror fesm-utils tsgen (TsGen.jl).

import * as Plot from "https://cdn.jsdelivr.net/npm/@observablehq/plot@0.6.16/+esm";
import katex from "https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/katex.mjs";

const root = document.getElementById("fb");
let D = null;                       // forcing.json
const S = { base: null, ov: new Map(), open: new Set(), runopts: "-r", outdir: "" };   // state: base name, overrides "g.k" -> value

// ---------------------------------------------------------------- resolution

const TEMPLATES = [[/^idx_\w+$/, "idx"], [/^snap_field_\w+$/, "snap_field"],
                   [/^snap_axis_\w+$/, "snap_axis"], [/^snap_\w+$/, "snap_snapshot"]];

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

function setOv(g, k, v) {
  const i = keyInfo(g, k);
  const id = g + "." + k;
  if (i && same(v, i.ref)) S.ov.delete(id); else S.ov.set(id, v);
  render();
}

// the timeline group ([ctrl] run_step)
const tgroup = () => { const r = eff("ctrl", "run_step", "ctrl"); return groupSource(r) ? r : "ctrl"; };
const climate = () => eff("comps", "climate", "snapshots");

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

// Known option lists, by "group.key" (or "template.key"); the idx method list
// comes from its comment.
const OPTIONS = {
  "ctrl.tstep_method": ["const", "rel", "cal"], "spinup.tstep_method": ["const", "rel", "cal"], "transient.tstep_method": ["const", "rel", "cal"],
  "comps.surface_method": ["chion", "smbpal", "smb_simple", "climate"],
  "sim.init_method": ["none", "equil", "recon", "recon_ref"],
  "barysealevel.method": ["const", "file", "fastiso"],
  "orbit.method": ["clock", "time", "fixed"], "ice_bnd.method": ["file", "scaled"],
  "snap.combine": ["anomaly", "absolute"], "snap_field.blend": ["linear", "ratio", "const"],
  "ocn_obm.mode": ["absolute", "anomaly"], "taper.outside": ["ref", "mean"],
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
// badge and the comment as help.
function field(g, k, opts = {}) {
  const info = keyInfo(g, k);
  if (!info) return h("div", { class: "fb-field fb-missing" }, `${g}.${k}: not in the schema`);
  const v = info.v, id = g + "." + k;
  let input;
  const options = opts.options || optionsOf(g, k, info);
  if (typeof v === "boolean") {
    input = h("input", { type: "checkbox", onchange: (e) => setOv(g, k, e.target.checked) });
    input.checked = v;
  } else if (options && !Array.isArray(v)) {
    input = h("select", { onchange: (e) => setOv(g, k, e.target.value) },
      ...options.map(o => h("option", { value: o }, o)), ...(options.includes(String(v)) ? [] : [h("option", { value: v }, String(v))]));
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
  const badge = h("span", { class: "fb-src fb-src-" + info.src, title: `value from the ${info.src}` + (info.src === "override" ? ` (${info.refsrc}: ${JSON.stringify(info.ref)})` : "") }, info.src);
  const help = (info.comment || "").split("\n")[0];
  return h("div", { class: "fb-field" },
    h("label", { title: info.comment, for: id }, h("code", {}, opts.label || k)),
    input, badge, help ? h("span", { class: "fb-help" }, help) : null);
}
function fields(g, keys, opts = {}) { return h("div", { class: "fb-fields" }, ...keys.map(k => field(g, k, opts[k] || {}))); }

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

function section(id, title, ...content) {
  return h("section", { class: "fb-section", id: "fb-" + id }, h("h3", {}, title), ...content);
}
function details(id, title, ...content) {
  const d = h("details", { class: "fb-details" }, h("summary", {}, title), ...content);
  if (S.open.has(id)) d.open = true;
  d.addEventListener("toggle", () => { if (d.open) S.open.add(id); else S.open.delete(id); });
  return d;
}
function note(text, cls = "fb-note") { return h("p", { class: cls }, text); }
function warn(text) { return note(text, "fb-warn"); }

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
function idxPar(g) {
  const p = {};
  for (const k of schemaKeys(g)) p[k] = eff(g, k);
  // normalisation is active only when both norm_t0 and norm_t1 are set in the group itself
  const def = groupDef(g);
  p.norm = !!(def && def.keys.norm_t0 && def.keys.norm_t1) || S.ov.has(g + ".norm_t0") || S.ov.has(g + ".norm_t1");
  return p;
}
function idxSeries(p) {
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
// f on the run's timeline (tl) for index group g; feedback methods: the envelope.
function indexOnRun(g, tl) {
  const p = idxPar(g), ser = idxSeries(p);
  const f = tl.tf.map(tf => applyWindow(p, tf, p.method === "series" ? (ser ? interp(ser, tf) : NaN) : evalOpen(p, tf - tl.tf0)));
  return { p, ser, f };
}

// ------------------------------------------------------------------- figures

const W = 640, H = 220;
function linePlot(series, { x = "t [yr]", y = "", rules = [], title = "", vrules = [], band = null, color = null }) {
  // series: [{label, t:[], f:[]}]
  const data = series.flatMap(s => s.t.map((t, i) => ({ t, f: s.f[i], label: s.label })));
  const marks = [Plot.ruleY([0], { stroke: "#999" })];
  if (band) marks.push(Plot.rectX([band], { x1: "x0", x2: "x1", fill: "#ffd", stroke: "none" }));
  marks.push(Plot.line(data, { x: "t", y: "f", stroke: "label", strokeWidth: 1.5 }));
  for (const r of rules) marks.push(Plot.ruleY([r.y], { stroke: r.stroke || "#c33", strokeDasharray: "4 3" }), Plot.text([{ t: data.length ? data[0].t : 0, f: r.y }], { x: "t", y: "f", text: () => r.label, dy: -6, textAnchor: "start", fill: r.stroke || "#c33" }));
  for (const r of vrules) marks.push(Plot.ruleX([r.x], { stroke: r.stroke || "#36c", strokeDasharray: "4 3" }), Plot.text([{ t: r.x }], { x: "t", text: () => r.label, frameAnchor: "top", dx: 4, textAnchor: "start", fill: r.stroke || "#36c" }));
  const ys = data.map(d => d.f).filter(Number.isFinite).concat(rules.map(r => r.y));
  let domain;
  if (ys.length) { const lo = Math.min(...ys), hi = Math.max(...ys); if (hi - lo < 1e-9) domain = [lo - 1, hi + 1]; }
  const el = Plot.plot({ width: W, height: H, marginLeft: 60, title, x: { label: x }, y: { label: y, grid: true, domain },
    color: series.length > 1 ? { legend: true } : undefined, marks });
  return h("figure", { class: "fb-fig" }, el);
}

// ------------------------------------------------------------------ sections

function secBase() {
  const sel = h("select", { onchange: (e) => { S.base = e.target.value; S.ov.clear(); render(); } },
    ...Object.keys(D.bases).map(b => h("option", { value: b }, b)));
  sel.value = S.base;
  const dom = eff("domain", "name", "?"), grid = eff("domain", "grid_hub", eff("domain", "grid_ice", ""));
  return section("base", "1. Base par file",
    h("div", { class: "fb-fields" }, h("div", { class: "fb-field" }, h("label", {}, h("code", {}, "-n")), sel,
      h("span", { class: "fb-help" }, `${base().file}: domain ${dom}${grid ? ", grid " + grid : ""}, climate backend ${climate()}`))),
    note("A base is one domain and climate backend; everything an experiment changes is a -p override of it. Changing the base clears the overrides."));
}

function secRun() {
  const tl = timeline();
  const T = tl.T, hasPhases = ["spinup", "transient"].some(g => groupSource(g) === "base");
  const stepOpts = ["ctrl", ...["spinup", "transient"].filter(g => groupSource(g) === "base")];
  const clockEq = tl.method === "const" ? `t_f = t_{\\mathrm{const}} = ${fmtK(tl.tc)}`
    : tl.method === "rel" ? `t_f = t` : `t_f = t_{\\mathrm{cal}} - ${fmtNum(tl.cref)}`;
  const content = [
    hasPhases ? fields("ctrl", ["run_step"], { run_step: { options: stepOpts } }) : null,
    fields(T, ["tstep_method", "tstep_const", "time_init", "time_end", "dtt"]),
    fields("ctrl", ["calendar", "calendar_ref", "restart"]),
    h("div", { class: "fb-eqrow" }, tex(`\\text{forcing clock:}\\quad ${clockEq}`),
      h("span", { class: "fb-help" }, "const: the model time advances, the forcing is held at tstep_const; rel: the forcing follows the model time (paleo records); cal: calendar years against calendar_ref")),
    linePlot([{ label: "t_f", t: tl.t, f: tl.tf }], { y: "forcing time t_f [yr]", title: "Forcing clock over the run" }),
  ];
  if (keyInfo(T, "time_ref")) content.push(h("h4", {}, "ESM periods"), fields(T, ["time_ref", "time_esm_ref", "time_hist", "time_proj", "clim_var", "clim_seed"]));
  content.push(h("h4", {}, "Components and output"),
    fields("comps", ["climate", "surface_method", "dt_clim", "with_ice_sheet", "with_isostasy", "with_climate", "with_surface", "with_shelf"], { climate: { options: [String(eff("comps", "climate"))] } }),
    note("The climate backend is fixed by the base: pick another base for another backend."),
    h("div", { class: "fb-fields" }, field("tm_1D", "dt", { label: "tm_1D.dt" }), field("tm_2D", "dt", { label: "tm_2D.dt" }), field("tm_rst", "dt", { label: "tm_rst.dt" })),
    h("h4", {}, "Initial state"),
    fields("sim", ["init_method", "relax", "opt"]));
  const im = String(eff("sim", "init_method", "none"));
  if (im.startsWith("recon")) content.push(fields("sim", ["recon_path", "recon_var", "recon_time", "recon_regions"]));
  if (im === "equil") content.push(fields("sim", ["init_equil_time"]));
  return section("run", "2. Run control", ...content);
}

const METHOD_EQ = {
  "const": (p) => `f = \\min\\!\\big(\\max(f_{\\min} + s\\,\\dot f_{\\max}\\,\\tau,\\ f_{\\min}),\\ f_{\\max}\\big)` ,
  "ramp-slope": (p) => `f = \\min\\!\\big(\\max(f_{\\mathrm{start}} + s\\,\\dot f_{\\max}\\,\\tau,\\ f_{\\min}),\\ f_{\\max}\\big)`,
  "ramp-time": (p) => `f = \\min\\!\\big(\\max(f_{\\mathrm{start}} + s\\,\\tfrac{|f_{\\max}-f_{\\min}|}{\\Delta t_{\\mathrm{ramp}}}\\,\\tau,\\ f_{\\min}),\\ f_{\\max}\\big)`,
  "ramp-time-step": (p) => `f = \\begin{cases} f_{\\mathrm{start}} + (f_{\\mathrm{ext}}-f_{\\mathrm{start}})\\,\\tau/\\Delta t_{\\mathrm{ramp}} & \\tau \\le \\Delta t_{\\mathrm{ramp}} \\\\ f_{\\mathrm{ext}} + (f_{\\mathrm{conv}}-f_{\\mathrm{ext}})\\,(\\tau-\\Delta t_{\\mathrm{ramp}})/\\Delta t_{\\mathrm{conv}} & \\tau \\le \\Delta t_{\\mathrm{ramp}}+\\Delta t_{\\mathrm{conv}} \\\\ f_{\\mathrm{conv}} & \\text{after} \\end{cases}`,
  "sin": (p) => `f = \\tfrac{f_{\\max}+f_{\\min}}{2} + \\tfrac{f_{\\max}-f_{\\min}}{2}\\,\\sin\\!\\big(2\\pi\\, t_e/\\Delta t_{\\mathrm{ramp}}\\big)`,
  "series": (p) => p.norm ? `f = \\frac{S(t_f) - S(t_0)}{S(t_1) - S(t_0)}` : `f = S(t_f)`,
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

function indexCard(name, tl) {
  const g = "idx_" + name;
  const src = groupSource(g);
  if (!src) return h("div", { class: "fb-card" }, warn(`[${g}]: no such group and no idx template`));
  const { p, ser, f } = indexOnRun(g, tl);
  const method = String(p.method);
  const isfb = FEEDBACK.includes(method);
  const keys = isfb ? FEEDBACK_KEYS : (METHOD_KEYS[method] || ["f_min", "f_max"]);
  const seriesFiles = Object.keys(D.series).sort();
  const opts = { series_file: { options: seriesFiles } };
  if (method === "series" && p.series_file && D.series[p.series_file] && D.series[p.series_file].kind === "netcdf")
    opts.series_var = { options: Object.keys(D.series[p.series_file].vars) };
  const content = [
    h("div", { class: "fb-cardhead" }, h("code", {}, `[${g}]`), h("span", { class: "fb-src fb-src-" + src }, src === "new" ? "new (template idx)" : src)),
    fields(g, ["method"]),
    fields(g, keys, opts),
    details("idx-extra-" + name, "noise, window, restarts", fields(g, ["sigma", "units", "restart_every_df", "win_t0", "win_t1", "win_taper", "win_fac", "win_ref", "win_above"])),
  ];
  // the equation with the values
  const s = num(p.df_sign) >= 0 ? 1 : -1;
  const vals = [];
  if (isfb) {
    content.push(tex(`\\text{${method}: response-driven, } \\dot f = \\text{controller}(\\dot V_{\\mathrm{ice}},\\ \\mathrm{tol}=${fmtNum(p.tol)}),\\quad |\\dot f| \\le ${fmtNum(p.df_dt_max)},\\quad f \\in [${fmtNum(p.f_min)},\\ ${fmtNum(p.f_max)}]`));
    content.push(note(`Feedback method: no closed form. The figure shows the fastest path, f_start + s·df_dt_max·τ, from f_start = ${s > 0 ? "f_min" : "f_max"} after dt_init = ${fmtNum(p.dt_init)} yr; the controller slows it while the ice volume responds faster than tol${p.with_kill ? ", and the run stops once the response equilibrates at a bound (with_kill)" : ""}.`));
  } else if (METHOD_EQ[method]) {
    content.push(tex(METHOD_EQ[method](p)));
    if (method === "series") {
      vals.push(`S: ${p.series_file}${p.series_var ? ":" + p.series_var : ""}`);
      if (p.norm) vals.push(`t_0 = ${asList(p.norm_t0).map(fmtNum).join("…")}, t_1 = ${asList(p.norm_t1).map(fmtNum).join("…")}` + (ser && ser.norm ? ` → S(t_0) = ${fmtNum(+ser.norm[0].toPrecision(4))}, S(t_1) = ${fmtNum(+ser.norm[1].toPrecision(4))}` : ""));
      if (!ser) vals.push("series not found in the exported data");
    } else {
      vals.push(`τ = max(0, t_f − t_f(t_init) − dt_init), dt_init = ${fmtNum(p.dt_init)} yr; s = ${s}; f_min = ${fmtNum(p.f_min)}, f_max = ${fmtNum(p.f_max)}`);
      if (method === "const" || method === "ramp-slope") vals.push(`df_dt_max = ${fmtNum(p.df_dt_max)} /yr` + (num(p.df_dt_max) > 0 ? ` → from ${fmtNum(s > 0 ? p.f_min : p.f_max)} to ${fmtNum(s > 0 ? p.f_max : p.f_min)} over ${fmtNum(Math.abs(num(p.f_max) - num(p.f_min)) / num(p.df_dt_max))} yr` : " → constant f = " + fmtNum(s > 0 ? p.f_min : p.f_max)));
      if (method === "ramp-time") vals.push(`dt_ramp = ${fmtNum(p.dt_ramp)} yr`);
      if (method === "ramp-time-step") vals.push(`dt_ramp = ${fmtNum(p.dt_ramp)} yr, dt_conv = ${fmtNum(p.dt_conv)} yr, f_conv = ${fmtNum(p.f_conv)}`);
      if (method === "sin") vals.push(`period dt_ramp = ${fmtNum(p.dt_ramp)} yr, t_e the time since the start of the run (held at 0 during dt_init)`);
      if (tl.method === "const") vals.push("the forcing clock is constant: τ = 0 and the index stays at f_start; use tstep_method = rel to advance it");
    }
  }
  if (num(p.win_t1) > num(p.win_t0)) vals.push(`window ${fmtNum(p.win_t0)} < t_f < ${fmtNum(p.win_t1)}: f → ${fmtNum(p.win_ref)} + ${fmtNum(p.win_fac)}·(f − ${fmtNum(p.win_ref)})${p.win_above ? " where f > win_ref" : ""}`);
  if (num(p.sigma) > 0) vals.push(`noise: + ε, ε ~ N(0, ${fmtNum(p.sigma)}²) every step`);
  if (vals.length) content.push(h("ul", { class: "fb-vals" }, ...vals.map(v => h("li", {}, v))));
  // figures
  const finite = f.filter(Number.isFinite);
  if (finite.length) {
    const rules = [];
    if (!isfb && method !== "series") rules.push({ y: num(p.f_min), label: "f_min" }, { y: num(p.f_max), label: "f_max" });
    content.push(linePlot([{ label: name, t: tl.t, f }], { y: `index ${name}`, rules, title: `Index ${name} over the run (model time)` }));
  }
  if (ser) {
    const vr = tl.method === "const" ? [{ x: tl.tc, label: "tstep_const" }] : [];
    if (p.norm) { for (const t of asList(p.norm_t0).map(num)) vr.push({ x: t, label: "t0", stroke: "#393" }); for (const t of asList(p.norm_t1).map(num)) vr.push({ x: t, label: "t1", stroke: "#393" }); }
    const band = tl.method !== "const" ? { x0: tl.tfmin, x1: tl.tfmax } : null;
    content.push(linePlot([{ label: "S", t: ser.time, f: ser.f }], { x: "forcing time t_f [yr]", y: p.norm ? "normalised series" : "series", vrules: vr, band, title: `The series (${ser.n} points${p.norm ? ", normalised" : ""})` + (ser.header ? ": " + ser.header.slice(0, 80) : "") }));
  }
  return h("div", { class: "fb-card" }, ...content);
}

function secIndices(tl) {
  const names = asList(eff("indices", "names", ""));
  const candidates = [...new Set([...Object.keys(base().groups), ...Object.keys(D.library.groups)])].filter(g => g.startsWith("idx_")).map(g => g.slice(4)).sort();
  return section("indices", "3. Indices",
    note("Named scalar series on the forcing clock, each a group [idx_<name>] from the par file, the forcing library or new from the idx template. Consumers (axes, modifiers, prescribed ice) refer to them by name."),
    namesEditor("indices", "names", candidates, "indices.names"),
    ...names.map(n => indexCard(n, tl)));
}

function secSnapshots(tl) {
  if (!["snapshots", "rembo"].includes(climate())) return null;
  const snaps = asList(eff("snap", "snapshots", ""));
  const axes = asList(eff("snap", "axes", ""));
  const fieldsL = asList(eff("snap", "fields", ""));
  const indices = asList(eff("indices", "names", ""));
  const candSnap = [...new Set([...Object.keys(base().groups), ...Object.keys(D.library.groups)])].filter(g => /^snap_\w+$/.test(g) && !/^snap_(axis|field)_/.test(g)).map(g => g.slice(5)).sort();
  const candAxes = [...new Set([...Object.keys(base().groups), ...Object.keys(D.library.groups)])].filter(g => /^snap_axis_/.test(g)).map(g => g.slice(10)).sort();
  const content = [
    note(climate() === "rembo" ? "REMBO supplies the atmosphere; the snapshot backend supplies the ocean from the reference snapshot." : "The climate is the reference snapshot plus, per field, the axes: each axis blends two snapshots (lo at 0, hi at 1) by the value of an index."),
    fields("snap", ["ref_name", "combine", "lapse", "f_p"]),
    namesEditor("snap", "snapshots", candSnap, "snap.snapshots"),
    h("div", { class: "fb-cards" }, ...snaps.map(s => {
      const g = "snap_" + s, src = groupSource(g);
      if (!src) return h("div", { class: "fb-card" }, warn(`[${g}]: no such group`));
      const bound = ["tas", "pr", "zs", "to", "so"].map(fld => { const v = asList(eff(g, fld, "")); return v.length ? h("li", {}, h("code", {}, fld), ": ", v.map(b => { const c = D.climdb[b]; return c ? `${b} (${c.filename.replace(/^fesmdata\//, "")}${Array.isArray(c.time) && c.time[0] > 1 ? ", " + c.time[0] + " slices" : ""})` : b; }).join("; ")) : null; });
      return h("div", { class: "fb-card" }, h("div", { class: "fb-cardhead" }, h("code", {}, `[${g}]`), h("span", { class: "fb-src fb-src-" + src }, src)),
        fields(g, ["exp", "monthly", "clim_sigma", "ocn_sigma"]), h("ul", { class: "fb-vals" }, ...bound));
    })),
    namesEditor("snap", "axes", candAxes, "snap.axes"),
    h("div", { class: "fb-cards" }, ...axes.map(a => {
      const g = "snap_axis_" + a, src = groupSource(g);
      if (!src) return h("div", { class: "fb-card" }, warn(`[${g}]: no such group`));
      const idx = String(eff(g, "index", ""));
      const w = [];
      if (idx && !indices.includes(idx)) w.push(warn(`index "${idx}" is not in indices.names`));
      for (const k of ["lo", "hi"]) { const v = String(eff(g, k, "")); if (v && !snaps.includes(v)) w.push(warn(`${k} = "${v}" is not in snap.snapshots`)); }
      return h("div", { class: "fb-card" }, h("div", { class: "fb-cardhead" }, h("code", {}, `[${g}]`), h("span", { class: "fb-src fb-src-" + src }, src)),
        fields(g, ["lo", "hi", "index"], { lo: { options: ["", ...snaps] }, hi: { options: snaps }, index: { options: indices } }), ...w);
    })),
    tex(`\\text{linear: } X = X_{\\mathrm{ref}} + \\sum_k a_k\\,(X_{\\mathrm{hi},k} - X_{\\mathrm{lo},k}),\\qquad \\text{ratio: } P = P_{\\mathrm{ref}} \\prod_k \\big(1 + a_k\\,(P_{\\mathrm{hi},k}/P_{\\mathrm{lo},k} - 1)\\big)`),
    note("a_k is the value of the axis index (lo = \"\": hi is itself the anomaly or ratio). Fields: " + fieldsL.map(f => `${f} (${eff("snap_field_" + f, "blend", "linear")}${asList(eff("snap_field_" + f, "axes", "")).length ? ", axes " + asList(eff("snap_field_" + f, "axes", "")).join(" ") : ""})`).join("; ") + ". Then the air temperature is brought from sea level to the surface, tas = tsl − z_srf·λ(m) with λ between lapse[1] (annual) and lapse[2] (summer), and the precipitation re-sensitised, pr = pr_sl · exp(β_p (tas − tsl)) with β_p = f_p."),
  ];
  // axes over the run
  const ax = axes.map(a => { const g = "snap_axis_" + a, idx = String(eff(g, "index", "")); return indices.includes(idx) ? { label: `${a} (index ${idx})`, t: tl.t, f: indexOnRun("idx_" + idx, tl).f } : null; }).filter(Boolean);
  if (ax.length) content.push(linePlot(ax, { y: "position on the axis", rules: [{ y: 0, label: "lo", stroke: "#393" }, { y: 1, label: "hi", stroke: "#393" }], title: "Axis positions over the run" }));
  return section("snapshots", "4. Snapshots and axes", ...content);
}

function secEsm(tl) {
  if (climate() !== "esm") return null;
  const T = tl.T;
  const periods = ["time_ref", "time_esm_ref", "time_hist", "time_proj"].map(k => { const v = asList(eff(T, k, "")).map(num); return v.length === 2 ? { k, x0: v[0], x1: v[1] } : null; }).filter(Boolean);
  if (tl.cal) periods.push({ k: "run (" + T + ")", x0: tl.method === "const" ? tl.tc : tl.t0, x1: tl.method === "const" ? tl.tc : tl.t1 });
  const fig = Plot.plot({ width: W, height: 160, marginLeft: 110, x: { label: "year CE", grid: true }, y: { label: null },
    marks: [Plot.barX(periods, { x1: "x0", x2: "x1", y: "k", fill: "#69b" }), Plot.ruleX([tl.method === "const" ? tl.tc : NaN], { stroke: "#c33", strokeDasharray: "4 3" })] });
  const useEsm = eff("esm", "use_esm", false), exp = String(eff("esm", "experiment", "ctrl"));
  const content = [
    note("The ESM backend applies a reference climatology at the current surface plus ESM anomalies relative to its reference period, following the periods of the timeline group (calendar years)."),
    fields("esm", ["experiment", "esm_name", "use_esm", "use_hist", "use_proj", "use_var", "basins", "par_file"]),
    fields("esm", ["lapse", "f_p", "f_ocn", "f_polar", "dT_threshold"]),
    h("figure", { class: "fb-fig" }, fig),
  ];
  if (useEsm && exp !== "ctrl") content.push(tex(`X(t) = X_{\\mathrm{ref}} + \\big(X_{\\mathrm{ESM}}(t) - \\overline{X_{\\mathrm{ESM}}}^{\\,\\mathrm{time\\_esm\\_ref}}\\big),\\quad \\text{historical to } ${fmtNum(asList(eff(T, "time_hist", [0, 0]))[1])},\\ \\text{projection from } ${fmtNum(asList(eff(T, "time_proj", [0, 0]))[0])}`),
    note("Air temperature and smb anomalies are differences, precipitation a ratio, the ocean annual means at the shelf base."));
  else if (!useEsm && exp === "1pctCO2") content.push(tex(`\\Delta T_a(t) = f_{\\mathrm{polar}}\\,\\min\\!\\big(17\\,\\mathrm{K}\\,\\tfrac{t - 1995}{305},\\ \\Delta T_{\\mathrm{thr}}\\big) = ${fmtNum(eff("esm", "f_polar"))}\\,\\min(\\ldots,\\ ${fmtNum(eff("esm", "dT_threshold"))}),\\qquad \\Delta T_o = f_{\\mathrm{ocn}}\\,\\Delta T_a = ${fmtNum(eff("esm", "f_ocn"))}\\,\\Delta T_a`),
    note("Homogeneous warming (use_esm = False with experiment 1pctCO2): a linear 17 K from 1995 to 2300, capped at dT_threshold, times f_polar."));
  else content.push(note(`experiment = "${exp}", use_esm = ${useEsm}: the reference climatology only (spin-up).`));
  return section("esm", "4. ESM forcing", ...content);
}

function secModifiers(tl) {
  const indices = asList(eff("indices", "names", ""));
  const cl = climate();
  const content = [note("Applied to the backend's climate in this fixed order; each an optional group with active.")];
  // anom
  const aActive = eff("anom", "active", false), aIdx = String(eff("anom", "index", ""));
  const anomEq = cl === "rembo"
    ? `\\Delta T_{\\mathrm{sum}} = f_{ta}\\,f(t) = ${fmtNum(eff("anom", "f_ta"))}\\,f,\\quad \\Delta T_{\\mathrm{ann}} = 1.3\\,\\Delta T_{\\mathrm{sum}},\\quad \\Delta T_o = f_{to}\\,\\Delta T_{\\mathrm{ann}} = ${fmtNum(eff("anom", "f_to"))}\\,\\Delta T_{\\mathrm{ann}}`
    : `\\Delta T_a = f_{ta}\\,f(t) = ${fmtNum(eff("anom", "f_ta"))}\\,f,\\quad \\Delta T_o = f_{to}\\,f = ${fmtNum(eff("anom", "f_to"))}\\,f,\\quad \\Delta S_o = f_{so}\\,f = ${fmtNum(eff("anom", "f_so"))}\\,f,\\quad P \\to P\\,e^{\\beta_p \\Delta T_a}`;
  const anomBody = [fields("anom", ["active", "index", "f_ta", "f_to", "f_so"], { index: { options: ["", ...indices] } })];
  if (aActive) {
    anomBody.push(tex(anomEq));
    if (!indices.includes(aIdx)) anomBody.push(warn(`anom.index = "${aIdx}" is not in indices.names`));
    else {
      const f = indexOnRun("idx_" + aIdx, tl).f;
      const ser = [];
      const fta = num(eff("anom", "f_ta")), fto = num(eff("anom", "f_to")), fso = num(eff("anom", "f_so"));
      if (cl === "rembo") { ser.push({ label: "dT_sum [K]", t: tl.t, f: f.map(x => x * fta) }, { label: "dT_ann [K]", t: tl.t, f: f.map(x => 1.3 * x * fta) }); if (fto) ser.push({ label: "dT_ocn [K]", t: tl.t, f: f.map(x => 1.3 * x * fta * fto) }); }
      else { ser.push({ label: "dTa [K]", t: tl.t, f: f.map(x => x * fta) }); if (fto) ser.push({ label: "dTo [K]", t: tl.t, f: f.map(x => x * fto) }); if (fso) ser.push({ label: "dSo [psu]", t: tl.t, f: f.map(x => x * fso) }); }
      anomBody.push(linePlot(ser, { y: "anomaly", title: `Homogeneous anomalies from index ${aIdx}` }));
    }
    if (cl === "esm" && String(eff("comps", "surface_method")) === "climate") anomBody.push(warn("anom is refused when the ESM supplies the smb (surface_method = climate)."));
  }
  content.push(h("div", { class: "fb-card" }, h("div", { class: "fb-cardhead" }, h("code", {}, "[anom]"), " homogeneous anomaly from an index"), ...anomBody));
  // ocn_fraction
  const oActive = eff("ocn_fraction", "active", false);
  const ocnBody = [fields("ocn_fraction", ["active", "f_to", "region"])];
  if (oActive) {
    ocnBody.push(tex(`T_o(z) = T_{o,\\mathrm{ref}}(z) + f_{to}\\,\\overline{\\big(T_{sl} - T_{sl,\\mathrm{ref}}\\big)}^{\\,\\mathrm{region}} = T_{o,\\mathrm{ref}} + ${fmtNum(eff("ocn_fraction", "f_to"))}\\,\\overline{\\Delta T_{sl}}`));
    if (aActive && num(eff("anom", "f_to")) !== 0) ocnBody.push(warn("ocn_fraction replaces the ocean temperature anomaly: set anom.f_to = 0."));
    if (cl !== "snapshots") ocnBody.push(warn("ocn_fraction needs ocean profiles and a sea-level air temperature: snapshots backend only."));
  }
  content.push(h("div", { class: "fb-card" }, h("div", { class: "fb-cardhead" }, h("code", {}, "[ocn_fraction]"), " ocean anomaly as a fraction of the air anomaly"), ...ocnBody));
  // taper, bias
  const tBody = [fields("taper", ["active", "weight", "var_defs", "outside"])];
  if (eff("taper", "active", false)) tBody.push(tex(`X = X_{\\mathrm{ref}} + w\\,(X - X_{\\mathrm{ref}})` + (String(eff("taper", "outside")) === "mean" ? ` + (1-w)\\,\\overline{(X - X_{\\mathrm{ref}})}^{\\,w=1}` : "")));
  content.push(h("div", { class: "fb-card" }, h("div", { class: "fb-cardhead" }, h("code", {}, "[taper]"), " the atmosphere towards the reference with a weight field"), ...tBody));
  const bBody = [fields("bias", ["active", "dtas", "fpr", "f_tas", "f_pr", "var_defs"])];
  if (eff("bias", "active", false)) bBody.push(tex(`T_a \\to T_a + ${fmtNum(eff("bias", "f_tas"))}\\,\\delta T_{\\mathrm{bias}},\\qquad P \\to P\\,\\big(1 + ${fmtNum(eff("bias", "f_pr"))}\\,(f_{P,\\mathrm{bias}} - 1)\\big)`));
  content.push(h("div", { class: "fb-card" }, h("div", { class: "fb-cardhead" }, h("code", {}, "[bias]"), " static bias correction of the atmosphere"), ...bBody));
  return section("modifiers", "5. Climate modifiers", ...content);
}

function secExtras(tl) {
  const indices = asList(eff("indices", "names", ""));
  const content = [];
  // sea level
  const m = String(eff("barysealevel", "method", "fastiso"));
  const slBody = [fields("barysealevel", ["method", "bsl_init", "sl_path", "sl_name"], { sl_path: { options: ["none", ...Object.keys(D.series).sort()] } })];
  if (m === "const") slBody.push(tex(`z_{sl} = ${fmtNum(eff("barysealevel", "bsl_init"))}\\ \\mathrm{m}`));
  else if (m === "fastiso") slBody.push(tex(`z_{sl}(t) = z_{sl,0} - \\Delta V_{\\mathrm{ice}}(t)/A_{\\mathrm{ocean}}`), note("Barystatic sea level from the ice volume of the FastIsostasy domains (bsl_init at the start)."));
  else if (m === "file") {
    const path = String(eff("barysealevel", "sl_path", "")), vname = String(eff("barysealevel", "sl_name", ""));
    const ser = seriesOf(path, vname);
    slBody.push(tex(`z_{sl}(t) = S_{sl}(t_f)`));
    if (ser) {
      slBody.push(linePlot([{ label: "z_sl [m]", t: tl.t, f: tl.tf.map(x => interp(ser, x)) }], { y: "sea level [m]", title: "Sea level over the run" }));
      slBody.push(linePlot([{ label: "S_sl", t: ser.time, f: ser.f }], { x: "forcing time t_f [yr]", y: "sea level [m]", title: `${path}: ${ser.n} points`, band: tl.method !== "const" ? { x0: tl.tfmin, x1: tl.tfmax } : null, vrules: tl.method === "const" ? [{ x: tl.tc, label: "tstep_const" }] : [] }));
    } else slBody.push(warn(`${path}${vname && vname !== "none" ? ":" + vname : ""}: not in the exported series (netCDF: set sl_name)`));
  }
  content.push(h("div", { class: "fb-card" }, h("div", { class: "fb-cardhead" }, h("code", {}, "[barysealevel]"), " sea level"), ...slBody));
  // orbit
  const om = String(eff("orbit", "method", "time"));
  const oBody = [fields("orbit", ["method", "time", "ecc", "obl", "per"])];
  oBody.push(om === "clock" ? tex(`\\text{orbit}(t_f): \\text{Laskar (2004) at the forcing clock}`)
    : om === "time" ? tex(`\\text{orbit}(${fmtK(num(eff("orbit", "time")))}\\ \\mathrm{yr}): \\text{fixed, Laskar (2004)}`)
    : tex(`e = ${fmtNum(eff("orbit", "ecc"))},\\ \\varepsilon = ${fmtNum(eff("orbit", "obl"))}^\\circ,\\ \\varpi = ${fmtNum(eff("orbit", "per"))}^\\circ`));
  content.push(h("div", { class: "fb-card" }, h("div", { class: "fb-cardhead" }, h("code", {}, "[orbit]"), " insolation of the surface models"), ...oBody));
  // ice_bnd
  if (groupSource("ice_bnd")) {
    const ib = [fields("ice_bnd", ["active", "method", "H_ice", "index", "slice_time", "var_defs"], { index: { options: ["", ...indices] } })];
    if (eff("ice_bnd", "active", false)) {
      const im = String(eff("ice_bnd", "method"));
      ib.push(im === "scaled" ? tex(`H_{\\mathrm{ref}} = f(t)\\,H_{\\mathrm{ice}}(${fmtK(num(eff("ice_bnd", "slice_time")))})`) : tex(`H_{\\mathrm{ref}} = H_{\\mathrm{ice}}(t_f)`));
      if (im === "scaled" && !indices.includes(String(eff("ice_bnd", "index")))) ib.push(warn("ice_bnd.index is not in indices.names"));
    }
    content.push(h("div", { class: "fb-card" }, h("div", { class: "fb-cardhead" }, h("code", {}, "[ice_bnd]"), " prescribed ice in the fixed and relaxation cells"), ...ib));
  }
  return section("extras", "6. Sea level, orbit, prescribed ice", ...content);
}

// --------------------------------------------------------------- the command

function fmtP(v) {
  if (typeof v === "boolean") return v ? "True" : "False";
  if (typeof v === "number") return fmtNum(v);
  const s = String(v);
  return /^[\w.\/{}:~&|,-]*$/.test(s) ? s : JSON.stringify(s);
}
function overrides() {
  // group order: the base's groups, then referenced groups
  const order = [...base().order];
  const add = (g) => { if (!order.includes(g)) order.push(g); };
  for (const n of asList(eff("indices", "names", ""))) add("idx_" + n);
  if (groupSource("snap")) { for (const n of asList(eff("snap", "snapshots", ""))) add("snap_" + n); for (const n of asList(eff("snap", "axes", ""))) add("snap_axis_" + n); }
  for (const id of S.ov.keys()) add(id.split(".")[0]);
  const out = [], blocks = [], problems = [];
  for (const g of order) {
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
  return { out, blocks, problems };
}
function fmtNml(v) {
  if (Array.isArray(v)) return v.map(fmtNml).join(" ");
  if (typeof v === "boolean") return v ? "True" : "False";
  if (typeof v === "number") return fmtNum(v);
  return JSON.stringify(String(v));
}
function parseP(text) {
  // runme -p tokens (whitespace separated; quoted values allowed) -> overrides
  const toks = text.match(/(?:[^\s"']+|"[^"]*"|'[^']*')+/g) || [];
  const errs = [];
  for (let tok of toks) {
    if (tok === "-p" || tok === "\\") continue;
    const m = tok.match(/^([\w]+)\.([\w]+)(==|=)(.*)$/s);
    if (!m) { errs.push(`not a key=value: ${tok}`); continue; }
    const [, g, k, op, raw] = m;
    const info = keyInfo(g, k);
    if (!info) { errs.push(`${g}.${k}: unknown group or key`); continue; }
    const items = (op === "==" ? raw.split(",") : [raw]).map(s => s.replace(/^["']|["']$/g, ""));
    const refel = Array.isArray(info.ref) ? (info.ref[0] ?? "") : info.ref;
    const conv = (s) => typeof refel === "boolean" ? /^(true|t|1|\.true\.)$/i.test(s) : typeof refel === "number" ? num(s) : s;
    const v = (op === "==" || Array.isArray(info.ref)) ? items.map(conv) : conv(items[0]);
    if (same(v, info.ref)) S.ov.delete(g + "." + k); else S.ov.set(g + "." + k, v);
  }
  return errs;
}

function secCommand() {
  const { out, blocks, problems } = overrides();
  const T = tgroup();
  const outdir = S.outdir || `output/${S.base.replace(/^yelmox_/, "").toLowerCase()}-run`;
  const cmd = [`runme ${S.runopts} -e yelmox`, `-n ${base().file}`, `-o ${outdir}`, ...(out.length ? ["-p " + out.join(" \\\n      ")] : [])].join(" \\\n   ");
  const pre = h("pre", { class: "fb-cmd" }, h("code", {}, cmd));
  const runopts = h("div", { class: "fb-fields" },
    h("div", { class: "fb-field" }, h("label", {}, h("code", {}, "runme options")), h("input", { type: "text", class: "fb-text", value: S.runopts, onchange: (e) => { S.runopts = e.target.value; render(); } }), h("span"), h("span", { class: "fb-help" }, "-r runs, -s stages, -rs submits; queue options as in the scripts (-q, -w, --omp)")),
    h("div", { class: "fb-field" }, h("label", {}, h("code", {}, "-o")), h("input", { type: "text", class: "fb-text", value: outdir, onchange: (e) => { S.outdir = e.target.value; render(); } }), h("span"), h("span", { class: "fb-help" }, "the run folder")));
  const copy = h("button", { class: "fb-btn", onclick: () => navigator.clipboard.writeText(cmd.replace(/ \\\n\s+/g, " ")) }, "copy (one line)");
  const ta = h("textarea", { class: "fb-ta", rows: 3, placeholder: "paste a -p list here, e.g. indices.names==at,ramp anom.active=True anom.index=ramp" });
  const load = h("button", { class: "fb-btn", onclick: () => { const errs = parseP(ta.value); render(); if (errs.length) alert(errs.join("\n")); } }, "load -p list");
  const reset = h("button", { class: "fb-btn", onclick: () => { S.ov.clear(); render(); } }, "reset to base");
  const content = [
    note("Only keys that differ from the base (or the library group or the schema default) are listed."),
    runopts, pre, h("div", {}, copy, " ", reset),
    ...problems.map(warn),
    ...(blocks.length ? [h("pre", { class: "fb-cmd" }, h("code", {}, blocks.join("\n\n")))] : []),
    h("h4", {}, "Load a -p list"), ta, h("div", {}, load),
  ];
  // the resolved configuration
  const order = [...base().order];
  for (const n of asList(eff("indices", "names", ""))) if (!order.includes("idx_" + n)) order.push("idx_" + n);
  if (groupSource("snap")) { for (const n of asList(eff("snap", "snapshots", ""))) if (!order.includes("snap_" + n)) order.push("snap_" + n); for (const n of asList(eff("snap", "axes", ""))) if (!order.includes("snap_axis_" + n)) order.push("snap_axis_" + n); }
  const rows = [];
  for (const g of order) {
    if (!groupSource(g)) continue;
    for (const k of schemaKeys(g)) {
      const i = keyInfo(g, k);
      rows.push(h("tr", { class: "fb-row-" + i.src }, h("td", {}, h("code", {}, `${g}.${k}`)), h("td", {}, h("code", {}, fmtNml(i.v))), h("td", {}, h("span", { class: "fb-src fb-src-" + i.src }, i.src)), h("td", { class: "fb-help" }, (i.comment || "").split("\n")[0])));
    }
  }
  content.push(details("resolved", `Resolved configuration (${rows.length} keys)`, h("table", { class: "fb-table" }, h("thead", {}, h("tr", {}, h("th", {}, "key"), h("th", {}, "value"), h("th", {}, "source"), h("th", {}, ""))), h("tbody", {}, ...rows))));
  return section("command", "7. The runme command", ...content);
}

// --------------------------------------------------------------------- render

function render() {
  const y = window.scrollY;
  root.replaceChildren();
  const tl = timeline();
  root.append(secBase(), secRun(), secIndices(tl), secSnapshots(tl), secEsm(tl), secModifiers(tl), secExtras(tl), secCommand());
  root.append(h("p", { class: "fb-help" }, `Data: forcing.json generated ${D.generated}${D.commit ? " from " + D.commit : ""} (scripts/forcing-builder/export.jl). Single-domain yelmox groups (no domain suffix); the snapshot fields of a run are not drawn.`));
  window.scrollTo(0, y);
}

fetch(new URL("forcing-builder/forcing.json", document.baseURI)).then(r => r.json()).then(d => {
  D = d;
  const q = new URLSearchParams(location.search);
  S.base = D.bases[q.get("base")] ? q.get("base") : (D.bases.yelmox_Greenland ? "yelmox_Greenland" : Object.keys(D.bases)[0]);
  if (q.get("p")) parseP(q.get("p"));
  render();
}).catch(err => { root.replaceChildren(h("p", { class: "fb-warn" }, "forcing.json could not be loaded: " + err)); });
