/* ===========================
 * Helpers: parsing & math
 * =========================== */
function parseMatrix(input) {
  const rows = input.trim().split(/\n|;+/).map(r => r.trim()).filter(Boolean);
  if (!rows.length) return [];
  const mat = rows.map(row => {
    const cols = row.split(/[,\s]+/).map(x => x.trim()).filter(Boolean).map(Number);
    if (cols.some(v => Number.isNaN(v))) throw new Error("Matrix contains non-numeric values.");
    return cols;
  });
  const w = mat[0].length;
  if (!mat.every(r => r.length === w)) throw new Error("Matrix rows have inconsistent lengths.");
  return mat;
}
function transpose(A){const m=A.length,n=A[0].length;const T=Array.from({length:n},()=>Array(m).fill(0));for(let i=0;i<m;i++)for(let j=0;j<n;j++)T[j][i]=A[i][j];return T;}
function matMul(A,B){const m=A.length,n=A[0].length,nB=B.length,k=B[0].length;if(n!==nB)throw new Error(`Dimension mismatch: (${m}x${n}) * (${nB}x${k}).`);const C=Array.from({length:m},()=>Array(k).fill(0));for(let i=0;i<m;i++){for(let j=0;j<k;j++){let s=0;for(let t=0;t<n;t++)s+=A[i][t]*B[t][j];C[i][j]=s;}}return C;}
function toHom(points){if(!points.length)return[];return points.map(p=>[...p,1]);}
function fromHom(pointsH){return pointsH.map(p=>{const w=p[p.length-1];const c=p.slice(0,-1);return Math.abs(w)>1e-12&&Math.abs(w-1)>1e-12?c.map(v=>v/w):c;});}
function applyTransform(points,T,useHom){
  if(!points.length) return [];
  const d=points[0].length;
  if(!points.every(p=>p.length===d)) throw new Error("All points must have same dimension.");
  if(useHom){
    if(!(T.length===d+1 && T[0].length===d+1)) throw new Error(`For homogeneous mode, transform must be ${(d+1)}x${(d+1)}.`);
    const R=transpose(matMul(T,transpose(toHom(points)))); return fromHom(R);
  }else{
    if(!(T.length===d && T[0].length===d)) throw new Error(`For linear mode, transform must be ${d}x${d}.`);
    return transpose(matMul(T,transpose(points)));
  }
}

/* ===========================
 * Presets
 * =========================== */
const dr = a => (a*Math.PI)/180;
function rot2D(theta){const t=dr(theta),c=Math.cos(t),s=Math.sin(t);return [[c,-s],[s,c]];}
function rot2DH(theta){const R=rot2D(theta);return [[R[0][0],R[0][1],0],[R[1][0],R[1][1],0],[0,0,1]];}
function rot3D(axis,theta){const t=dr(theta),c=Math.cos(t),s=Math.sin(t); if(axis==='x')return [[1,0,0],[0,c,-s],[0,s,c]]; if(axis==='y')return [[c,0,s],[0,1,0],[-s,0,c]]; return [[c,-s,0],[s,c,0],[0,0,1]];}
function rot3DH(axis,theta){const R=rot3D(axis,theta);return [[R[0][0],R[0][1],R[0][2],0],[R[1][0],R[1][1],R[1][2],0],[R[2][0],R[2][1],R[2][2],0],[0,0,0,1]];}
function scaleH2(sx,sy){return [[sx,0,0],[0,sy,0],[0,0,1]];}
function scaleH3(sx,sy,sz){return [[sx,0,0,0],[0,sy,0,0],[0,0,sz,0],[0,0,0,1]];}
function translateH2(tx,ty){return [[1,0,tx],[0,1,ty],[0,0,1]];}
function translateH3(tx,ty,tz){return [[1,0,0,tx],[0,1,0,ty],[0,0,1,tz],[0,0,0,1]];}
function shear2DH(shx,shy){return [[1,shx,0],[shy,1,0],[0,0,1]];}

/* ===========================
 * DOM refs
 * =========================== */
const el = {
  btn2d: document.getElementById('btn-2d'),
  btn3d: document.getElementById('btn-3d'),
  pointsLabel: document.getElementById('points-label'),
  transformLabel: document.getElementById('transform-label'),
  pointsText: document.getElementById('points-text'),
  transformText: document.getElementById('transform-text'),
  resetPoints: document.getElementById('reset-points'),
  resetTransform: document.getElementById('reset-transform'),
  chkHom: document.getElementById('chk-hom'),
  angle: document.getElementById('angle'),
  angleVal: document.getElementById('angle-val'),
  sx: document.getElementById('sx'), sy: document.getElementById('sy'), sz: document.getElementById('sz'),
  tx: document.getElementById('tx'), ty: document.getElementById('ty'), tz: document.getElementById('tz'),
  shx: document.getElementById('shx'), shy: document.getElementById('shy'),
  btnRotation: document.getElementById('btn-rotation'),
  btnScale: document.getElementById('btn-scale'),
  btnTranslation: document.getElementById('btn-translation'),
  btnShear: document.getElementById('btn-shear'),
  axisWrap: document.getElementById('axis-wrap'),
  axisSelect: document.getElementById('axis-select'),
  szWrap: document.getElementById('sz-wrap'),
  tzWrap: document.getElementById('tz-wrap'),
  shear2d: document.getElementById('shear-2d'),
  polyOpts: document.getElementById('poly-opts'),
  chkConnect: document.getElementById('chk-connect'),
  chkClose: document.getElementById('chk-close'),
  error: document.getElementById('error'),
  svg2d: document.getElementById('svg2d'),
  three3d: document.getElementById('three3d'),
  outOrig: document.getElementById('out-original'),
  outTrans: document.getElementById('out-transformed'),
};

let mode = '2D'; // '2D' | '3D'

/* ===========================
 * Defaults
 * =========================== */
const defaults = {
  points2D: `-1 -1
1 -1
1 1
-1 1`,
  transform2D_H: `1 0 1.2
0 1 0.4
0 0 1`,
  transform2D_L: `1 0
0 1`,
  points3D: `-1 -1 -1
1 -1 -1
1 1 -1
-1 1 -1
-1 -1 1
1 -1 1
1 1 1
-1 1 1`,
  transform3D_H: `1 0 0 0.8
0 1 0 0.4
0 0 1 0.2
0 0 0 1`,
  transform3D_L: `1 0 0
0 1 0
0 0 1`,
};

function setMode(newMode){
  mode = newMode;
  el.btn2d.classList.toggle('selected', mode==='2D');
  el.btn3d.classList.toggle('selected', mode==='3D');

  const useHom = el.chkHom.checked;
  el.pointsLabel.textContent = mode==='2D' ? 'Points (N×2)' : 'Points (N×3)';
  el.transformLabel.textContent = mode==='2D'
    ? (useHom?'Transform (3×3)':'Transform (2×2)')
    : (useHom?'Transform (4×4)':'Transform (3×3)');

  // show/hide 2D/3D-only controls
  el.axisWrap.classList.toggle('hidden', mode==='2D');
  el.szWrap.classList.toggle('hidden', mode==='2D');
  el.tzWrap.classList.toggle('hidden', mode==='2D');
  el.shear2d.classList.toggle('hidden', mode!=='2D');
  el.polyOpts.classList.toggle('hidden', mode!=='2D');

  el.svg2d.classList.toggle('hidden', mode!=='2D');
  el.three3d.classList.toggle('hidden', mode!=='3D');

  // if user switches modes, ensure we have matching example defaults
  if(mode==='2D' && !el.pointsText.value.trim()) el.pointsText.value = defaults.points2D;
  if(mode==='3D' && !el.pointsText.value.trim()) el.pointsText.value = defaults.points3D;
  if(mode==='2D') el.transformText.value = useHom ? defaults.transform2D_H : defaults.transform2D_L;
  else el.transformText.value = useHom ? defaults.transform3D_H : defaults.transform3D_L;

  renderAll();
}

function updateLabels(){
  const useHom = el.chkHom.checked;
  el.transformLabel.textContent = mode==='2D'
    ? (useHom?'Transform (3×3)':'Transform (2×2)')
    : (useHom?'Transform (4×4)':'Transform (3×3)');
}

/* ===========================
 * 2D Rendering (SVG)
 * =========================== */
function render2D(original, transformed){
  const svg = el.svg2d;
  const w = 560, h = 420, pad = 24;

  // clear
  while(svg.firstChild) svg.removeChild(svg.firstChild);

  // bounds
  const all = original.concat(transformed);
  const xs = all.map(p=>p[0]), ys = all.map(p=>p[1]);
  const minX = Math.min(-5, ...xs), maxX = Math.max(5, ...xs);
  const minY = Math.min(-5, ...ys), maxY = Math.max(5, ...ys);
  const spanX = (maxX-minX)||10, spanY = (maxY-minY)||10;
  const sx = x => pad + ((x-minX)/spanX)*(w-2*pad);
  const sy = y => h - (pad + ((y-minY)/spanY)*(h-2*pad));

  const line = (x1,y1,x2,y2,cls) => {
    const L = document.createElementNS("http://www.w3.org/2000/svg","line");
    L.setAttribute("x1",x1);L.setAttribute("y1",y1);L.setAttribute("x2",x2);L.setAttribute("y2",y2);
    L.setAttribute("class",cls);
    svg.appendChild(L);
  };

  // grid
  for(let i=0;i<=10;i++){
    const gx=minX+(i/10)*spanX, gy=minY+(i/10)*spanY;
    line(sx(gx), sy(minY), sx(gx), sy(maxY), "grid");
    line(sx(minX), sy(gy), sx(maxX), sy(gy), "grid");
  }
  // axes
  line(sx(0), sy(minY), sx(0), sy(maxY), "axis");
  line(sx(minX), sy(0), sx(maxX), sy(0), "axis");

  const connect = el.chkConnect.checked;
  const closePoly = el.chkClose.checked;

  const pathStr = pts => pts.map((p,i)=>`${i?'L':'M'}${sx(p[0]).toFixed(2)},${sy(p[1]).toFixed(2)}`).join(' ') + (closePoly && pts.length>2 ? ' Z':'');

  // original
  if(connect && original.length>1){
    const P=document.createElementNS("http://www.w3.org/2000/svg","path");
    P.setAttribute("d",pathStr(original)); P.setAttribute("class","orig");
    svg.appendChild(P);
  }
  original.forEach(p=>{
    const C=document.createElementNS("http://www.w3.org/2000/svg","circle");
    C.setAttribute("cx",sx(p[0])); C.setAttribute("cy",sy(p[1])); C.setAttribute("r","4");
    C.setAttribute("class","orig");
    svg.appendChild(C);
  });

  // transformed
  if(connect && transformed.length>1){
    const P=document.createElementNS("http://www.w3.org/2000/svg","path");
    P.setAttribute("d",pathStr(transformed)); P.setAttribute("class","trans");
    svg.appendChild(P);
  }
  transformed.forEach(p=>{
    const R=document.createElementNS("http://www.w3.org/2000/svg","rect");
    R.setAttribute("x",sx(p[0])-4); R.setAttribute("y",sy(p[1])-4); R.setAttribute("width","8"); R.setAttribute("height","8");
    R.setAttribute("class","trans");
    svg.appendChild(R);
  });
}

/* ===========================
 * 3D Rendering (Three.js)
 * =========================== */
let threeCtx = null; // { scene, camera, renderer, controls, groupOrig, groupTrans }
function ensureThree(){
  if (!window.THREE || !window.THREE.OrbitControls) return null;
  if (threeCtx) return threeCtx;

  const container = el.three3d;
  // clear
  container.innerHTML = '';
  const renderer = new THREE.WebGLRenderer({ antialias:true });
  renderer.setSize(container.clientWidth, container.clientHeight);
  container.appendChild(renderer.domElement);

  const scene = new THREE.Scene();
  const camera = new THREE.PerspectiveCamera(55, container.clientWidth/container.clientHeight, 0.1, 100);
  camera.position.set(3.5, 2.8, 3.5);

  const controls = new THREE.OrbitControls(camera, renderer.domElement);

  scene.add(new THREE.GridHelper(10,10, 0x9ca3af, 0xe5e7eb));
  scene.add(new THREE.AxesHelper(3));

  const amb = new THREE.AmbientLight(0xffffff, 0.8);
  const dir = new THREE.DirectionalLight(0xffffff, 0.9);
  dir.position.set(5,5,5);
  scene.add(amb); scene.add(dir);

  const groupOrig = new THREE.Group();
  const groupTrans = new THREE.Group();
  scene.add(groupOrig); scene.add(groupTrans);

  function animate(){
    requestAnimationFrame(animate);
    groupOrig.rotation.y += 0.0025;
    groupTrans.rotation.y += 0.0025;
    renderer.render(scene, camera);
  }
  animate();

  threeCtx = { scene, camera, renderer, controls, groupOrig, groupTrans, container };
  window.addEventListener('resize', ()=>{
    if (!threeCtx) return;
    const { container, camera, renderer } = threeCtx;
    const w = container.clientWidth, h = container.clientHeight || 1;
    camera.aspect = w / h; camera.updateProjectionMatrix();
    renderer.setSize(w, h);
  });
  return threeCtx;
}

function render3D(original, transformed){
  const ctx = ensureThree();
  const container = el.three3d;

  if (!ctx) {
    container.innerHTML = '<div class="center subtle" style="display:grid;place-items:center;height:100%;">3D preview needs Three.js.</div>';
    return;
  }
  const { groupOrig, groupTrans } = ctx;
  // clear groups
  while(groupOrig.children.length) groupOrig.remove(groupOrig.children[0]);
  while(groupTrans.children.length) groupTrans.remove(groupTrans.children[0]);

  const mkSphere = (color=0x2563eb) => {
    const geo = new THREE.SphereGeometry(0.07, 16, 16);
    const mat = new THREE.MeshStandardMaterial({ color });
    return new THREE.Mesh(geo, mat);
  };

  original.forEach(p=>{
    const m=mkSphere(0x2563eb); m.position.set(p[0]||0,p[1]||0,p[2]||0); groupOrig.add(m);
  });
  transformed.forEach(p=>{
    const m=mkSphere(0xdc2626); m.position.set(p[0]||0,p[1]||0,p[2]||0); groupTrans.add(m);
  });
}

/* ===========================
 * Glue
 * =========================== */
function currentPointsDefault(){ return mode==='2D' ? defaults.points2D : defaults.points3D; }
function currentTransformDefault(){
  const hom = el.chkHom.checked;
  if (mode==='2D') return hom ? defaults.transform2D_H : defaults.transform2D_L;
  return hom ? defaults.transform3D_H : defaults.transform3D_L;
}

function pasteMatrix(mat){
  el.transformText.value = mat.map(r=>r.join(' ')).join('\n');
  renderAll();
}

function compute(){
  try {
    el.error.classList.add('hidden');
    const P = parseMatrix(el.pointsText.value || currentPointsDefault());
    const T = parseMatrix(el.transformText.value || currentTransformDefault());
    const Q = applyTransform(P, T, el.chkHom.checked);
    return { P, Q };
  } catch(e){
    el.error.textContent = e.message || String(e);
    el.error.classList.remove('hidden');
    return { P: [], Q: [] };
  }
}

function renderAll(){
  const { P, Q } = compute();

  // numeric
  el.outOrig.textContent = P.map(r=>r.join('\t')).join('\n') || '—';
  el.outTrans.textContent = Q.map(r=>r.join('\t')).join('\n') || '—';

  // draw
  if (mode==='2D') render2D(P, Q);
  else render3D(P, Q);
}

/* ===========================
 * Events
 * =========================== */
el.btn2d.addEventListener('click', ()=> setMode('2D'));
el.btn3d.addEventListener('click', ()=> setMode('3D'));

el.chkHom.addEventListener('change', ()=>{
  updateLabels();
  // swap example transform to correct size when toggled
  el.transformText.value = currentTransformDefault();
  renderAll();
});

el.resetPoints.addEventListener('click', ()=>{
  el.pointsText.value = currentPointsDefault();
  renderAll();
});
el.resetTransform.addEventListener('click', ()=>{
  el.transformText.value = currentTransformDefault();
  renderAll();
});

['points-text','transform-text'].forEach(id=>{
  document.getElementById(id).addEventListener('input', renderAll);
});

el.chkConnect.addEventListener('change', renderAll);
el.chkClose.addEventListener('change', renderAll);

// presets
el.angle.addEventListener('input', ()=> el.angleVal.textContent = el.angle.value);
el.btnRotation.addEventListener('click', ()=>{
  const angle = parseFloat(el.angle.value);
  if (mode==='2D') pasteMatrix(el.chkHom.checked ? rot2DH(angle) : rot2D(angle));
  else {
    const axis = el.axisSelect.value;
    pasteMatrix(el.chkHom.checked ? rot3DH(axis, angle) : rot3D(axis, angle));
  }
});
el.btnScale.addEventListener('click', ()=>{
  const sx = parseFloat(el.sx.value), sy=parseFloat(el.sy.value);
  if (mode==='2D') pasteMatrix(el.chkHom.checked ? scaleH2(sx,sy) : [[sx,0],[0,sy]]);
  else {
    const sz=parseFloat(el.sz.value);
    pasteMatrix(el.chkHom.checked ? scaleH3(sx,sy,sz) : [[sx,0,0],[0,sy,0],[0,0,sz]]);
  }
});
el.btnTranslation.addEventListener('click', ()=>{
  const tx=parseFloat(el.tx.value), ty=parseFloat(el.ty.value);
  if (mode==='2D') pasteMatrix(el.chkHom.checked ? translateH2(tx,ty) : [[1,0],[0,1]]);
  else {
    const tz=parseFloat(el.tz.value);
    pasteMatrix(el.chkHom.checked ? translateH3(tx,ty,tz) : [[1,0,0],[0,1,0],[0,0,1]]);
  }
});
el.btnShear.addEventListener('click', ()=>{
  const shx=parseFloat(el.shx.value), shy=parseFloat(el.shy.value);
  pasteMatrix(el.chkHom.checked ? shear2DH(shx,shy) : [[1,shx],[shy,1]]);
});

// initial values & render
el.pointsText.value = defaults.points2D;
el.transformText.value = defaults.transform2D_H;
setMode('2D');
