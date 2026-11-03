"""
claims_lint: a mechanical check that what a repository SAYS is backed by what it CONTAINS.
Generalises the CNC project's doc_review.py (forbidden claims + explicit exemptions) and adds two checks that catch the failures found in the 2026-11 audit:

  1. TEXT RULES      forbidden / review-needed claims in docs and code comments (EN + HU), negation- and quotation-aware, exemption only with a written reason
  2. CODE RULES      (Python, AST) a verdict word (PROVED / CERTIFIED / VERIFIED ...) printed or returned UNCONDITIONALLY in a function;
                     a metric (cpu_load, power, latency, safety_score, accuracy ...) assigned a numeric literal; a rigged random draw next to 'leak/violation'
  3. LEDGER          claims_ledger.json: every listed claim must have (a) an evidence file that exists, (b) each number of the claim present in the evidence file,
                     (c) its documentation snippets present in the named docs and containing the same numbers  -> a doc/log mismatch (e.g. 60,803 vs 100,803) fails

Severity:  FAIL = blocks release;  REVIEW = needs a human decision recorded as an exemption;  INFO = exempted (shown, never hidden).
Exemption (same or previous line):   claims:allow(RULE_ID) reason of at least 15 characters         (a missing/short reason is itself a FAIL)
Quoted old claim block (e.g. CHANGES.md):  a line containing  claims:quote  exempts that line and the following contiguous block.
Usage:  python claims_lint.py <repo_dir> [--json out.json] [--ledger claims_ledger.json] [--no-code]
Exit code 0 only if there is no FAIL and no REVIEW.
"""
from __future__ import annotations
import ast, json, re, sys
from dataclasses import dataclass, asdict
from pathlib import Path

TEXT_EXT = {".md", ".html", ".htm", ".tex", ".txt", ".rst", ".py", ".js", ".ts", ".tsx", ".bio", ".yml", ".yaml", ".json", ".toml"}
SKIP_DIRS = {".git", "node_modules", "venv", ".venv", "site-packages", "__pycache__", "target", "build", "dist", ".next", "vendor", "_claims_fixtures"}
MAX_BYTES = 1_500_000
ALLOW = re.compile(r"claims:allow\(([A-Z0-9\-]+)\)\s*(.*)")
QUOTE = "claims:quote"
NEGATION = re.compile(r"\b(not|no|never|without|nor|neither|nem|nincs|soha|sem|n't|yet|until|withdrawn|visszavon\w*|non)\b", re.I)
QUOTED_CONTEXT = re.compile(r"(earlier|previous|old|former|prior)\s+(version|claim|text|statement|readme)|withdrawn|korábbi|régi|visszavon|was\s+(wrong|unsupported)|were\s+(wrong|unsupported|not\s+supported)", re.I)


@dataclass(frozen=True)
class Rule:
    id: str
    severity: str          # FAIL | REVIEW
    pattern: str
    why: str
    instead: str

    @property
    def rx(self):
        return re.compile(self.pattern, re.I)


RULES = [
    Rule("CERT-CLAIM", "FAIL",
         r"(SIL[-\s]?[1-4]\s*(certified|compliant|tanúsít\w*)|TÜV[-\s]?(certified|tanúsít\w*)|certified\s+(product|device|safe\w*|stage\s*\d?|by)|DO-178C\s+certified|FDA[-\s]?(cleared|approved)|FDA\s+510\(k\)\s+(cleared|approved)|tanúsított\s+(termék|rendszer|szoftver)|(NATO|DoD)[-\s/]+(ready|compliant|certified)|ready\s+for\s+(NATO|DoD|SIL[-\s]?[1-4])\b)",
         "States a certification or regulatory status that no external body has granted.",
         "State the target and the real status: 'roadmap toward ...; not assessed'."),
    Rule("ABSOLUTE", "FAIL",
         r"(100\s?%\s*(?:[\w-]+\s+){0,2}(safe\w*|secur\w*|proof|proven|certain\w*|accura\w*|coverage|protect\w*|safety|determinis\w*)|mathematical(ly)?\s+(certain\w*|impossib\w*|guarantee\w*)|matematikai(lag)?\s+(bizonyoss\w*|lehetetlen\w*|garant\w*)|unhackable|un-hackable|bulletproof|golyóálló|unbreakable|cannot\s+crash|can(no|')t\s+be\s+(hacked|bypassed|defeated)|immune\s+to|zero\s+(false[-\s]?positive\w*|jitter)|nulla\s+téves|hacking\s+is\s+a\s+mathematical\s+impossibility)",
         "Absolute claim about safety, security or accuracy that a finite test or proof cannot support.",
         "Name the property, the model and the evidence ('proved for the float32 step function, see ...')."),
    Rule("PERF-NUMBER", "FAIL",
         r"(\b0\.0005\s*ms\b|\b5\s*[–-]\s*20\s*ns\b|\b1600\s*[x×]|<\s*10\s*ns\b|\b800\s*mW\b|\b0\.5\s*mW\b|\bnanosecond\w*\s+(response|latency|guarantee)|nanoszekundum\w*\s+(válasz|késleltet|garant))",
         "A specific performance figure that recurs in earlier materials and was found to come from scripted constants.",
         "Use a measured value with its log, or delete. Cite the ledger id via claims:allow(PERF-NUMBER)."),
    Rule("MEDICAL-REG", "REVIEW",
         r"(\bFDA\b|510\(k\)|\bPMA\b|IEC\s*62304|Class\s+(II|III|C)\s+(device|medical|software)|patient\s+safety|implant\w*)",
         "Regulatory or medical context: any claim here needs regulatory sign-off.",
         "Record the sign-off as an exemption, or remove the claim."),
    Rule("PROOF-OVERSTATE", "REVIEW",
         r"(mathematically\s+proven|formally\s+verified\s+(system|device|logic|binary)|provably\s+(safe|secure|correct)|matematikailag\s+bizonyított|proved\s+correct\s+for\s+all)",
         "A proof claim: needs the model, the property and the solver evidence next to it.",
         "Say what was proved, about which model, with which tool (see the ledger)."),
]

VERDICT_WORDS = re.compile(r"(PROVED|CERTIFIED|VERIFIED|Mathematical Certainty|MATHEMATICALLY PROVEN|SOUND|BULLETPROOF)", re.I)
METRIC_NAME = re.compile(r"^(cpu_load\w*|power\w*|energy\w*|latency\w*|safety_score\w*|accuracy\w*|false_positive\w*|violations_leaked\w*|jitter\w*)$", re.I)


@dataclass
class Finding:
    file: str
    line: int
    rule: str
    severity: str
    text: str
    note: str = ""


def iter_files(root: Path):
    for p in sorted(root.rglob("*")):
        if p.is_dir() or any(part in SKIP_DIRS for part in p.parts):
            continue
        if p.suffix.lower() in TEXT_EXT and p.stat().st_size <= MAX_BYTES and not p.name.endswith((".min.js", ".lock")) and p.name not in ("claims_ledger.json", "claims_lint.py", "test_claims_lint.py"):  # the linter's own rule text and tests are not documents
            yield p


def exemptions(lines):
    """line numbers (1-based) -> (rule_id, reason) from allow markers, and the set of quoted lines"""
    allow, bad, quoted = {}, [], set()
    active = seen = False
    for i, ln in enumerate(lines, 1):
        m = ALLOW.search(ln)
        if m:
            rid, reason = m.group(1), m.group(2).strip(" -->*/#")
            if len(reason) < 15: bad.append((i, rid))
            for t in (i, i + 1): allow.setdefault(t, set()).add(rid)
        if QUOTE in ln: active, seen = True, False
        if active:
            quoted.add(i)
            if ln.strip(): seen = True
            elif seen: active = False
    return allow, bad, quoted


def text_findings(path: Path, rel: str):
    out = []
    try: text = path.read_text(encoding="utf-8", errors="replace")
    except OSError: return out
    lines = text.splitlines()
    allow, bad, quoted = exemptions(lines)
    for ln, rid in bad:
        out.append(Finding(rel, ln, "EXEMPTION-REASON", "FAIL", lines[ln - 1].strip()[:140], f"claims:allow({rid}) needs a written reason of at least 15 characters"))
    for rule in RULES:
        for m in rule.rx.finditer(text):
            ln = text.count("\n", 0, m.start()) + 1
            line = lines[ln - 1] if ln <= len(lines) else ""
            sev = rule.severity
            note = rule.why
            if rule.id in allow.get(ln, ()):
                sev, note = "INFO", "exempted by claims:allow"
            elif ln in quoted or QUOTED_CONTEXT.search(line):
                sev, note = "INFO", "quoted/withdrawn context"
            elif NEGATION.search(re.split(r'[.;!?]\s', text[max(0, m.start() - 90):m.start()])[-1]):
                sev, note = "INFO", "negated nearby"
            out.append(Finding(rel, ln, rule.id, sev, line.strip()[:160], note))
    return out


class CodeVisitor(ast.NodeVisitor):
    def __init__(self, rel, lines):
        self.rel, self.lines, self.f, self.stack = rel, lines, [], []

    def _flag(self, node, rule, sev, note):
        ln = getattr(node, "lineno", 0)
        self.f.append(Finding(self.rel, ln, rule, sev, self.lines[ln - 1].strip()[:140] if 0 < ln <= len(self.lines) else "", note))

    def generic_visit(self, node):
        self.stack.append(node)
        super().generic_visit(node)
        self.stack.pop()

    def _conditional(self):
        for n in reversed(self.stack):
            if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef, ast.Module)): return False
            if isinstance(n, (ast.If, ast.IfExp, ast.While, ast.For, ast.ExceptHandler, ast.Try, ast.Match)): return True
        return False

    def visit_Call(self, node):
        is_print = isinstance(node.func, ast.Name) and node.func.id == "print"
        if is_print:
            consts = [c.value for a in node.args for c in ast.walk(a) if isinstance(c, ast.Constant) and isinstance(c.value, str)]
            if any(VERDICT_WORDS.search(s) for s in consts) and not self._conditional():
                self._flag(node, "HARD-VERDICT", "REVIEW", "a verdict word is printed unconditionally: is it computed here or just asserted?")
        self.generic_visit(node)

    def visit_Return(self, node):
        if isinstance(node.value, ast.Dict) or isinstance(node.value, ast.Constant):
            consts = [c.value for c in ast.walk(node.value) if isinstance(c, ast.Constant) and isinstance(c.value, str)]
            if any(VERDICT_WORDS.search(s) for s in consts) and not self._conditional():
                self._flag(node, "HARD-VERDICT", "REVIEW", "a verdict word is returned unconditionally")
        self.generic_visit(node)

    def visit_Assign(self, node):
        if isinstance(node.value, ast.Constant) and isinstance(node.value.value, (int, float)) and not isinstance(node.value.value, bool):
            for t in node.targets:
                name = t.attr if isinstance(t, ast.Attribute) else getattr(t, "id", "")
                if METRIC_NAME.match(name or "") and node.value.value not in (0, 0.0):
                    self._flag(node, "METRIC-CONSTANT", "REVIEW", f"'{name}' set to the literal {node.value.value}: measurement or assumption? (label it as an assumption, or compute it)")
        self.generic_visit(node)

    def visit_If(self, node):
        uses_random = any(isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute) and n.func.attr in ("random", "uniform", "rand") for n in ast.walk(node.test))
        if uses_random:
            end = getattr(node, "end_lineno", node.lineno)
            body = " ".join(self.lines[node.lineno - 1:end])
            if re.search(r"leak|violation|miss|fail|breach", body, re.I):
                self._flag(node, "RIGGED-RANDOM", "REVIEW", "a random draw decides a safety outcome: the result comes from the script, not from the system")
        self.generic_visit(node)


def code_findings(path: Path, rel: str):
    try: src = path.read_text(encoding="utf-8", errors="replace"); tree = ast.parse(src)
    except (SyntaxError, OSError, ValueError): return []
    allow, _, quoted = exemptions(src.splitlines())
    v = CodeVisitor(rel, src.splitlines()); v.visit(tree)
    for f in v.f:
        if f.rule in allow.get(f.line, ()): f.severity, f.note = "INFO", "exempted by claims:allow"
    return v.f


def norm(s: str) -> str:
    return re.sub(r"(?<=\d)[,  ](?=\d{3}\b)", "", s)


def ledger_findings(root: Path, ledger_path: Path):
    out = []
    try: data = json.loads(ledger_path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as e:
        return [Finding(str(ledger_path), 0, "LEDGER", "FAIL", "", f"cannot read ledger: {e}")]
    statuses = {"PROVED", "MEASURED", "ESTIMATED", "SIMULATED", "ROADMAP", "NOT_DONE"}
    for c in data.get("claims", []):
        cid = c.get("id", "?"); rel = str(ledger_path.relative_to(root)) if ledger_path.is_relative_to(root) else str(ledger_path)
        if c.get("status") not in statuses: out.append(Finding(rel, 0, "LEDGER", "FAIL", cid, f"status must be one of {sorted(statuses)}"))
        ev_text = ""
        for ev in c.get("evidence", []):
            p = root / ev["file"]
            if not p.is_file(): out.append(Finding(rel, 0, "LEDGER-EVIDENCE", "FAIL", cid, f"evidence file missing: {ev['file']}")); continue
            t = p.read_text(encoding="utf-8", errors="replace"); ev_text += norm(t)
            if ev.get("contains") and not re.search(ev["contains"], t):
                out.append(Finding(rel, 0, "LEDGER-EVIDENCE", "FAIL", cid, f"pattern not found in {ev['file']}: {ev['contains']}"))
        if c.get("status") in ("PROVED", "MEASURED", "ESTIMATED", "SIMULATED") and not c.get("evidence"):
            out.append(Finding(rel, 0, "LEDGER-EVIDENCE", "FAIL", cid, "a PROVED/MEASURED/ESTIMATED/SIMULATED claim needs at least one evidence file"))
        for n in c.get("numbers", []):
            ev_tok = norm(str(n["evidence"] if isinstance(n, dict) else n)); doc_tok = norm(str(n["doc"] if isinstance(n, dict) else n))
            if ev_tok not in ev_text: out.append(Finding(rel, 0, "LEDGER-NUMBER", "FAIL", cid, f"number {ev_tok} does not occur in the evidence files"))
            if isinstance(n, dict):
                try:
                    dec = len(doc_tok.split(".")[1]) if "." in doc_tok else 0
                    if abs(round(float(ev_tok), dec) - float(doc_tok)) > 10 ** (-dec) * 0.51:
                        out.append(Finding(rel, 0, "LEDGER-NUMBER", "FAIL", cid, f"doc value {doc_tok} is not a rounding of the evidence value {ev_tok}"))
                except ValueError:
                    out.append(Finding(rel, 0, "LEDGER-NUMBER", "FAIL", cid, f"cannot compare {ev_tok} with {doc_tok}"))
        for d in c.get("docs", []):
            p = root / d["file"]
            if not p.is_file(): out.append(Finding(rel, 0, "LEDGER-DOC", "FAIL", cid, f"doc missing: {d['file']}")); continue
            t = p.read_text(encoding="utf-8", errors="replace")
            if d["text"] not in t:
                out.append(Finding(rel, 0, "LEDGER-DOC", "FAIL", cid, f"snippet not found in {d['file']}: {d['text'][:60]}")); continue
            for n in c.get("numbers", []):
                doc_tok = norm(str(n["doc"] if isinstance(n, dict) else n))
                if doc_tok not in norm(d["text"]): out.append(Finding(rel, 0, "LEDGER-NUMBER", "FAIL", cid, f"number {doc_tok} missing from the doc snippet in {d['file']}"))
    return out


def run(root: Path, ledger: Path | None = None, code: bool = True):
    findings = []
    for p in iter_files(root):
        rel = str(p.relative_to(root))
        findings += text_findings(p, rel)
        if code and p.suffix == ".py": findings += code_findings(p, rel)
    lp = ledger or (root / "claims_ledger.json")
    if lp.is_file(): findings += ledger_findings(root, lp)
    return findings


def main(argv):
    if not argv: print(__doc__); return 2
    root = Path(argv[0]).resolve(); ledger = None; jp = None; code = True
    if "--ledger" in argv: ledger = Path(argv[argv.index("--ledger") + 1])
    if "--json" in argv: jp = Path(argv[argv.index("--json") + 1])
    if "--no-code" in argv: code = False
    fs = run(root, ledger, code)
    cnt = {s: sum(1 for f in fs if f.severity == s) for s in ("FAIL", "REVIEW", "INFO")}
    for sev in ("FAIL", "REVIEW"):
        for f in [x for x in fs if x.severity == sev]:
            print(f"{sev:6s} {f.rule:16s} {f.file}:{f.line}  {f.text}\n       -> {f.note}")
    print(f"\n{root.name}: FAIL {cnt['FAIL']} | REVIEW {cnt['REVIEW']} | INFO (exempted) {cnt['INFO']}")
    if jp: jp.write_text(json.dumps([asdict(f) for f in fs], ensure_ascii=False, indent=1), encoding="utf-8")
    return 1 if cnt["FAIL"] or cnt["REVIEW"] else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
