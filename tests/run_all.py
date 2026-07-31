"""Run every TACO test file and return a non-zero exit code if any failed.

The individual files are self-executing scripts rather than a pytest suite, so
there was no single command that covered all of them and nothing would have
caught a regression in one file while another passed.

Run from the repo root:  python tests/run_all.py
No external tools required.
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
FILES = [
    "test_concordance.py",
    "test_purify.py",
    "test_review_fixes.py",
    "test_v135_rescue.py",
]

if __name__ == "__main__":
    failed = []
    counts = {}
    for name in FILES:
        path = os.path.join(HERE, name)
        if not os.path.isfile(path):
            print("MISSING  %s" % name)
            failed.append(name)
            continue
        res = subprocess.run([sys.executable, path], capture_output=True,
                             text=True)
        n_pass = res.stdout.count("PASS  ")
        counts[name] = n_pass
        if res.returncode != 0:
            failed.append(name)
            print("FAILED   %-24s (%d passed)" % (name, n_pass))
            # Only the failures are worth the reader's screen space.
            for line in res.stdout.splitlines():
                if line.startswith(("FAIL", "ERROR")):
                    print("    " + line)
            if res.stderr.strip():
                print("    stderr: " + res.stderr.strip().splitlines()[-1])
        else:
            print("ok       %-24s (%d passed)" % (name, n_pass))

    total = sum(counts.values())
    print("\n%d tests across %d files" % (total, len(counts)))
    print("ALL TESTS PASSED" if not failed
          else "%d FILE(S) FAILED: %s" % (len(failed), ", ".join(failed)))
    sys.exit(1 if failed else 0)
