import assert from "node:assert/strict"
import { execFileSync } from "node:child_process"
import {
  mkdtempSync,
  rmSync,
  unlinkSync,
  writeFileSync,
} from "node:fs"
import os from "node:os"
import path from "node:path"
import test from "node:test"

import {
  sourceFiles,
  treeDigest,
} from "../scripts/generate_frontend_notices.mjs"


test("source hashing ignores tracked files deleted from the working tree", () => {
  const root = mkdtempSync(path.join(os.tmpdir(), "openbio-notices-"))

  try {
    execFileSync("git", ["init", "--quiet"], { cwd: root })
    writeFileSync(path.join(root, "deleted.ts"), "old source\n")
    execFileSync("git", ["-c", "core.autocrlf=false", "add", "deleted.ts"], {
      cwd: root,
    })
    unlinkSync(path.join(root, "deleted.ts"))
    writeFileSync(path.join(root, "current.ts"), "current source\n")

    const files = sourceFiles(root)

    assert.deepEqual(files, ["current.ts"])
    assert.match(treeDigest(root, files), /^[a-f0-9]{64}$/)
  } finally {
    rmSync(root, { recursive: true, force: true })
  }
})
