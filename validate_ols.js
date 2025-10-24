// Lightweight OLS validator for the provided CSV using Node.js
// Model: z = b0 + b1*x + b2*y with
//   z = pMIC(Mtb-H37Rv), x = pKD-NTD, y = pKD-FL

const fs = require('fs');

function parseNumber(s) {
  if (s == null) return NaN;
  s = String(s).trim();
  if (!s) return NaN;
  // Treat qualifiers like '>10' or '<0.1' as missing
  if (s.startsWith('>') || s.startsWith('<')) return NaN;
  // Remove thousands separators if any (not expected here)
  const cleaned = s.replace(/,/g, '');
  const v = Number(cleaned);
  return Number.isFinite(v) ? v : NaN;
}

function readCSV(path) {
  const raw = fs.readFileSync(path, 'utf8');
  // Handle CR, LF, or CRLF line endings
  const lines = raw.split(/\r\n|\n|\r/).filter((l) => l.trim().length > 0);
  if (lines.length === 0) throw new Error('Empty CSV');
  const header = lines[0].split(',');
  const rows = lines.slice(1).map((line) => line.split(','));
  return { header, rows };
}

function selectColumns(header, rows, cols) {
  const idx = cols.map((name) => {
    const i = header.indexOf(name);
    if (i === -1) throw new Error(`Column not found: ${name}`);
    return i;
  });
  const data = rows.map((r) => idx.map((i) => (i < r.length ? r[i] : '')));
  return data;
}

function buildDesignMatrix(data) {
  // data rows: [z, x, y]
  const xs = [];
  const zs = [];
  for (const [zRaw, xRaw, yRaw] of data) {
    const z = parseNumber(zRaw);
    const x = parseNumber(xRaw);
    const y = parseNumber(yRaw);
    if (Number.isFinite(z) && Number.isFinite(x) && Number.isFinite(y)) {
      xs.push([1, x, y]);
      zs.push(z);
    }
  }
  return { X: xs, z: zs };
}

function transpose(A) {
  const m = A.length, n = A[0].length;
  const AT = Array.from({ length: n }, () => Array(m).fill(0));
  for (let i = 0; i < m; i++) {
    for (let j = 0; j < n; j++) AT[j][i] = A[i][j];
  }
  return AT;
}

function matMul(A, B) {
  const m = A.length, n = A[0].length, p = B[0].length;
  const C = Array.from({ length: m }, () => Array(p).fill(0));
  for (let i = 0; i < m; i++) {
    for (let k = 0; k < n; k++) {
      const aik = A[i][k];
      for (let j = 0; j < p; j++) C[i][j] += aik * B[k][j];
    }
  }
  return C;
}

function matVec(A, v) {
  const m = A.length, n = A[0].length;
  const out = Array(m).fill(0);
  for (let i = 0; i < m; i++) {
    let s = 0;
    for (let j = 0; j < n; j++) s += A[i][j] * v[j];
    out[i] = s;
  }
  return out;
}

function invert3(M) {
  // Inverse of 3x3 matrix using adjugate/determinant
  const [a, b, c] = M[0];
  const [d, e, f] = M[1];
  const [g, h, i] = M[2];
  const A = e * i - f * h;
  const B = -(d * i - f * g);
  const C = d * h - e * g;
  const D = -(b * i - c * h);
  const E = a * i - c * g;
  const F = -(a * h - b * g);
  const G = b * f - c * e;
  const H = -(a * f - c * d);
  const I = a * e - b * d;
  const det = a * A + b * B + c * C;
  if (Math.abs(det) < 1e-12) throw new Error('Singular matrix');
  const inv = [
    [A / det, D / det, G / det],
    [B / det, E / det, H / det],
    [C / det, F / det, I / det],
  ];
  return inv;
}

function mean(v) {
  return v.reduce((a, b) => a + b, 0) / v.length;
}

function linAlgSolve(X, z) {
  const XT = transpose(X);
  const XTX = matMul(XT, X);
  const XTZ = matVec(XT, z);
  const XTXinv = invert3(XTX);
  const beta = matVec(XTXinv, XTZ);
  return { beta, XTXinv };
}

function computeMetrics(X, z, beta, XTXinv) {
  const n = X.length;
  const k = 2; // number of predictors
  const yhat = X.map((row) => row[0] * beta[0] + row[1] * beta[1] + row[2] * beta[2]);
  const resid = z.map((zi, i) => zi - yhat[i]);
  const zbar = mean(z);
  const ssRes = resid.reduce((s, r) => s + r * r, 0);
  const ssTot = z.reduce((s, zi) => s + (zi - zbar) ** 2, 0);
  const r2 = 1 - ssRes / ssTot;
  const dfRes = n - (k + 1);
  const mseRes = ssRes / dfRes;
  const mseModel = (ssTot - ssRes) / k;
  const fStat = mseRes > 0 ? mseModel / mseRes : Infinity;
  const logLike = -0.5 * n * (Math.log(2 * Math.PI) + Math.log(mseRes) + 1);
  const aic = -2 * logLike + 2 * (k + 1);
  const bic = -2 * logLike + Math.log(n) * (k + 1);
  const se = [0, 0, 0].map((_, j) => Math.sqrt(mseRes * XTXinv[j][j]));
  return { yhat, resid, r2, ssRes, ssTot, dfRes, mseRes, fStat, logLike, aic, bic, se };
}

function main() {
  const path = process.argv[2] || 'StarDrop-for-3D-least-squares-data-LEIs-AtlasCBS-columns-added-test-1.csv';
  const { header, rows } = readCSV(path);
  const zName = 'pMIC(Mtb-H37Rv)';
  const xName = 'pKD-NTD';
  const yName = 'pKD-FL';
  const data = selectColumns(header, rows, [zName, xName, yName]);
  const { X, z } = buildDesignMatrix(data);
  const n = X.length;
  if (n < 3) throw new Error('Not enough valid rows after filtering');
  const { beta, XTXinv } = linAlgSolve(X, z);
  const { r2, dfRes, fStat, logLike, aic, bic, se } = computeMetrics(X, z, beta, XTXinv);

  console.log('OLS validation (z = b0 + b1*x + b2*y)');
  console.log(`Rows used: ${n}`);
  console.log(`Columns: z='${zName}', x='${xName}', y='${yName}'`);
  console.log('Coefficients:');
  console.log(`  const (b0): ${beta[0].toFixed(4)}  (SE ${se[0].toFixed(3)})`);
  console.log(`  b1   (x):  ${beta[1].toFixed(4)}  (SE ${se[1].toFixed(3)})`);
  console.log(`  b2   (y):  ${beta[2].toFixed(4)}  (SE ${se[2].toFixed(3)})`);
  console.log('Model fit:');
  console.log(`  R^2: ${r2.toFixed(3)}  Adj R^2: ${(1 - (1 - r2) * ((n - 1) / dfRes)).toFixed(3)}`);
  console.log(`  F-statistic: ${fStat.toFixed(3)}`);
  console.log(`  Log-Likelihood: ${logLike.toFixed(3)}  AIC: ${aic.toFixed(3)}  BIC: ${bic.toFixed(3)}`);
}

if (require.main === module) {
  try {
    main();
  } catch (e) {
    console.error('Error:', e.message);
    process.exit(1);
  }
}
