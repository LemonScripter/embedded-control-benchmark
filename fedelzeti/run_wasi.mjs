// Futtatja a wasm32-wasip1 binarist Node WASI-val: node run_wasi.mjs <wasm> <tv-fajl> <iteraciok>
import { readFileSync } from 'node:fs';
import { WASI } from 'node:wasi';
import { argv } from 'node:process';
const [, , wasmPath, tv, iters] = argv;
const wasi = new WASI({ version: 'preview1', args: ['kolumban_fedelzeti', tv, iters ?? '20000'], preopens: { '.': '.' }, returnOnExit: true });
const mod = await WebAssembly.compile(readFileSync(wasmPath));
const inst = await WebAssembly.instantiate(mod, wasi.getImportObject());
wasi.start(inst);
