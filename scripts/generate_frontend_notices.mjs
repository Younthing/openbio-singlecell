#!/usr/bin/env node

import { execFileSync } from "node:child_process"
import { createHash } from "node:crypto"
import {
  copyFileSync,
  existsSync,
  readFileSync,
  readdirSync,
  statSync,
  writeFileSync,
} from "node:fs"
import path from "node:path"
import { fileURLToPath } from "node:url"

const scriptDir = path.dirname(fileURLToPath(import.meta.url))
const pluginRoot = path.dirname(scriptDir)

function usage() {
  console.log(`Usage: node scripts/generate_frontend_notices.mjs [options]

Generate the license files that must accompany the paired OpenBio frontend dist.

Options:
  --frontend-root <path>  ComfyUI_frontend checkout (default: sibling checkout)
  --dist <path>           Built frontend directory (default: <frontend>/dist)
  --check                 Verify existing files without writing
  -h, --help              Show this help

Set OPENBIO_PNPM_CLI to pnpm.cjs when pnpm is not available through corepack.`)
}

function parseArgs(argv) {
  let frontendRoot = path.resolve(pluginRoot, "..", "ComfyUI_frontend")
  let dist
  let check = false

  for (let i = 0; i < argv.length; i += 1) {
    const arg = argv[i]
    if (arg === "-h" || arg === "--help") {
      usage()
      process.exit(0)
    } else if (arg === "--check") {
      check = true
    } else if (arg === "--frontend-root" || arg === "--dist") {
      const value = argv[++i]
      if (!value) throw new Error(`${arg} requires a path`)
      if (arg === "--frontend-root") frontendRoot = path.resolve(value)
      else dist = path.resolve(value)
    } else {
      throw new Error(`Unknown argument: ${arg}`)
    }
  }

  return { frontendRoot, dist: dist ?? path.join(frontendRoot, "dist"), check }
}

function normalizeText(value) {
  return value.replace(/^\uFEFF/, "").replace(/\r\n/g, "\n").trim()
}

function sha256(value) {
  return createHash("sha256").update(value).digest("hex")
}

function formatPerson(value) {
  if (!value) return ""
  if (typeof value === "string") return value
  return [value.name, value.email, value.url].filter(Boolean).join("; ")
}

function sourceUrl(pkg) {
  const repository =
    typeof pkg.repository === "string" ? pkg.repository : pkg.repository?.url
  const source = String(pkg.homepage ?? repository ?? "")
    .replace(/^git\+/, "")
    .replace(/\.git$/, "")
  if (/^[\w.-]+\/[\w.-]+$/.test(source)) return `https://github.com/${source}`
  return source.replace(/^git@github\.com:/, "https://github.com/")
}

function escapeCell(value) {
  return String(value || "—")
    .replace(/\|/g, "\\|")
    .replace(/\r?\n/g, " ")
}

function compareText(a, b) {
  return a < b ? -1 : a > b ? 1 : 0
}

function licenseFiles(packageRoot) {
  return readdirSync(packageRoot)
    .filter((name) =>
      /^(licen[sc]e|copying|copyright|notice)([._-].*)?$/i.test(name),
    )
    .sort(compareText)
}

function runPnpmLicenses(frontendRoot) {
  const args = [
    "--filter",
    "@comfyorg/comfyui-frontend...",
    "licenses",
    "list",
    "--prod",
    "--json",
  ]
  const explicitCli = process.env.OPENBIO_PNPM_CLI
  const activeCli = process.env.npm_execpath?.includes("pnpm")
    ? process.env.npm_execpath
    : undefined
  let stdout

  if (explicitCli || activeCli) {
    stdout = execFileSync(process.execPath, [explicitCli ?? activeCli, ...args], {
      cwd: frontendRoot,
      encoding: "utf8",
      maxBuffer: 64 * 1024 * 1024,
    })
  } else {
    if (process.platform === "win32") {
      stdout = execFileSync(
        process.env.ComSpec ?? "cmd.exe",
        ["/d", "/s", "/c", ["corepack", "pnpm", ...args].join(" ")],
        {
          cwd: frontendRoot,
          encoding: "utf8",
          maxBuffer: 64 * 1024 * 1024,
        },
      )
    } else {
      stdout = execFileSync("corepack", ["pnpm", ...args], {
        cwd: frontendRoot,
        encoding: "utf8",
        maxBuffer: 64 * 1024 * 1024,
      })
    }
  }

  return JSON.parse(stdout.slice(stdout.indexOf("{")))
}

export function sourceFiles(frontendRoot) {
  try {
    return execFileSync(
      "git",
      ["-C", frontendRoot, "ls-files", "--cached", "--others", "--exclude-standard", "-z"],
    )
      .toString("utf8")
      .split("\0")
      .filter(Boolean)
      .filter((relative) => {
        const candidate = path.join(frontendRoot, relative)
        return existsSync(candidate) && statSync(candidate).isFile()
      })
  } catch {
    return walkFiles(frontendRoot, new Set([".git", "dist", "node_modules"]))
  }
}

function walkFiles(root, excludedNames = new Set(), current = "") {
  const files = []
  for (const entry of readdirSync(path.join(root, current), { withFileTypes: true })) {
    if (excludedNames.has(entry.name)) continue
    const relative = path.join(current, entry.name)
    if (entry.isDirectory()) files.push(...walkFiles(root, excludedNames, relative))
    else if (entry.isFile()) files.push(relative)
  }
  return files
}

export function treeDigest(root, files) {
  const hash = createHash("sha256")
  for (const relative of [...files].sort(compareText)) {
    hash.update(relative.replace(/\\/g, "/"))
    hash.update("\0")
    hash.update(readFileSync(path.join(root, relative)))
    hash.update("\0")
  }
  return hash.digest("hex")
}

function baselineCommit() {
  const manifest = JSON.parse(readFileSync(path.join(pluginRoot, "release_manifest.json"), "utf8"))
  return manifest.comfyui_frontend.commit
}

function detectUnknownLicense(texts) {
  const joined = texts.join("\n")
  if (
    joined.includes("Permission is hereby granted, free of charge") &&
    joined.includes('THE SOFTWARE IS PROVIDED "AS IS"')
  ) {
    return "MIT"
  }
  if (joined.includes("Apache License") && joined.includes("Version 2.0")) {
    return "Apache-2.0"
  }
  if (joined.includes("Permission to use, copy, modify, and/or distribute this software")) {
    return "ISC"
  }
  return "Unknown (license text supplied by package)"
}

function collectPackages(grouped) {
  const packages = new Map()

  for (const [license, items] of Object.entries(grouped)) {
    for (const item of items) {
      for (const packageRoot of item.paths) {
        const pkg = JSON.parse(
          readFileSync(path.join(packageRoot, "package.json"), "utf8"),
        )
        const key = `${pkg.name}@${pkg.version}`
        if (packages.has(key)) continue

        const texts = licenseFiles(packageRoot).map((filename) => ({
          filename,
          text: normalizeText(readFileSync(path.join(packageRoot, filename), "utf8")),
        }))
        packages.set(key, {
          key,
          name: pkg.name,
          version: pkg.version,
          license: license === "Unknown" ? detectUnknownLicense(texts.map((x) => x.text)) : license,
          author: formatPerson(pkg.author),
          source: sourceUrl(pkg),
          texts,
        })
      }
    }
  }

  return [...packages.values()].sort((a, b) => compareText(a.key, b.key))
}

function addFallbackTexts(packages, frontendLicense) {
  const allTexts = packages.flatMap((pkg) => pkg.texts.map((entry) => entry.text))
  const mitSource = allTexts.find((text) =>
    text.includes("Permission is hereby granted, free of charge"),
  )
  const apacheSource = allTexts.find(
    (text) => text.includes("Apache License") && text.includes("Version 2.0"),
  )
  const iscSource = allTexts.find((text) =>
    text.includes("Permission to use, copy, modify, and/or distribute this software"),
  )
  if (!mitSource || !apacheSource || !iscSource) {
    throw new Error("Could not locate canonical MIT, Apache-2.0, and ISC terms")
  }

  const mitTerms = mitSource.slice(
    mitSource.indexOf("Permission is hereby granted, free of charge"),
  )
  const fallbacks = {
    MIT: `MIT License\n\nThe installed package archive omitted its license file. Package identity and upstream source are recorded in the inventory.\n\n${mitTerms}`,
    "Apache-2.0": apacheSource,
    "GPL-3.0-only": frontendLicense,
    ISC: iscSource,
  }
  const missing = []

  for (const pkg of packages) {
    if (pkg.texts.length) continue
    const text = fallbacks[pkg.license]
    if (!text) {
      missing.push(`${pkg.key} (${pkg.license})`)
      continue
    }
    pkg.usesCanonicalTerms = true
    pkg.texts.push({ filename: `canonical ${pkg.license} terms`, text })
  }

  if (missing.length) {
    throw new Error(`Packages without distributable license terms:\n${missing.join("\n")}`)
  }
}

function fenceFor(text) {
  const longest = Math.max(0, ...[...text.matchAll(/`+/g)].map((match) => match[0].length))
  return "`".repeat(Math.max(3, longest + 1))
}

function renderNotice({ frontendRoot, dist, packages, frontendLicense, upstreamNotice }) {
  const packageJson = JSON.parse(readFileSync(path.join(frontendRoot, "package.json"), "utf8"))
  const lockHash = sha256(readFileSync(path.join(frontendRoot, "pnpm-lock.yaml")))
  const sourceHash = treeDigest(frontendRoot, sourceFiles(frontendRoot))
  const distFiles = walkFiles(dist).filter(
    (relative) => relative !== "LICENSE" && relative !== "THIRD_PARTY_NOTICES.md",
  )
  const distHash = treeDigest(dist, distFiles)
  const canonicalTermsCount = packages.filter((pkg) => pkg.usesCanonicalTerms).length
  const texts = new Map()
  const licenseCounts = new Map()

  for (const pkg of packages) {
    licenseCounts.set(pkg.license, (licenseCounts.get(pkg.license) ?? 0) + 1)
    pkg.textIds = pkg.texts.map(({ filename, text }) => {
      const id = `license-${sha256(text).slice(0, 16)}`
      const current = texts.get(id)
      if (current && current.text !== text) throw new Error(`License hash collision: ${id}`)
      if (!current) texts.set(id, { id, text, sources: [] })
      texts.get(id).sources.push(`${pkg.key} (${filename})`)
      return id
    })
  }

  const lines = [
    "# ComfyUI_frontend third-party notices",
    "",
    "This file accompanies the prebuilt OpenBio frontend. It inventories the installed production dependency graph and reproduces the license and notice files available in those packages. Exact duplicate texts are stored once and referenced by ID.",
    "",
    `- Frontend: @comfyorg/comfyui-frontend ${packageJson.version}`,
    `- Source baseline commit: ${baselineCommit()}`,
    `- Source tree SHA-256: ${sourceHash}`,
    `- Frontend dist SHA-256 (excluding this notice and adjacent LICENSE): ${distHash}`,
    `- pnpm lockfile SHA-256: ${lockHash}`,
    `- Production package/version pairs: ${packages.length}`,
    `- Unique license/notice texts: ${texts.size}`,
    `- Packages using canonical SPDX terms because their installed archive omitted a license file: ${canonicalTermsCount}`,
    "",
    "The frontend's own GPL-3.0-only terms are in the adjacent `LICENSE` file. Internal workspace packages built from the same source tree are covered by that source distribution; this inventory covers installed production packages reported by pnpm. When an installed package omits a license file, its available author/source metadata remains in the inventory and canonical terms for its declared SPDX license are included below. This cannot reconstruct a copyright line omitted by the upstream archive; release reviewers must consult the linked source when additional notices are required.",
    "",
    "## License summary",
    "",
    "| License | Packages |",
    "| --- | ---: |",
    ...[...licenseCounts.entries()]
      .sort(([a], [b]) => compareText(a, b))
      .map(([license, count]) => `| ${escapeCell(license)} | ${count} |`),
    "",
    "## Package inventory",
    "",
    "| Package | License | Author | Source | Text IDs |",
    "| --- | --- | --- | --- | --- |",
    ...packages.map(
      (pkg) =>
        `| ${escapeCell(pkg.key)} | ${escapeCell(pkg.license)} | ${escapeCell(pkg.author)} | ${escapeCell(pkg.source)} | ${pkg.textIds.map((id) => `[${id}](#${id})`).join(", ")} |`,
    ),
    "",
    "## Upstream project notices",
    "",
    upstreamNotice || "No separate upstream notice text was present.",
    "",
    "## License and notice texts",
    "",
  ]

  for (const { id, text, sources } of [...texts.values()].sort((a, b) =>
    compareText(a.id, b.id),
  )) {
    const fence = fenceFor(text)
    lines.push(
      `### ${id}`,
      "",
      `Referenced by ${sources.length} package file${sources.length === 1 ? "" : "s"}.`,
      "",
      `${fence}text`,
      text,
      fence,
      "",
    )
  }

  return `${lines.join("\n").trimEnd()}\n`
}

function main() {
  const { frontendRoot, dist, check } = parseArgs(process.argv.slice(2))
  const required = [
    path.join(frontendRoot, "package.json"),
    path.join(frontendRoot, "pnpm-lock.yaml"),
    path.join(frontendRoot, "node_modules"),
    path.join(frontendRoot, "LICENSE"),
    path.join(dist, "index.html"),
  ]
  for (const item of required) {
    if (!existsSync(item)) throw new Error(`Required path is missing: ${item}`)
  }

  const frontendLicense = normalizeText(readFileSync(path.join(frontendRoot, "LICENSE"), "utf8"))
  const upstreamNoticePath = path.join(frontendRoot, "THIRD_PARTY_NOTICES.md")
  const upstreamNotice = existsSync(upstreamNoticePath)
    ? normalizeText(readFileSync(upstreamNoticePath, "utf8"))
    : ""
  const packages = collectPackages(runPnpmLicenses(frontendRoot))
  addFallbackTexts(packages, frontendLicense)
  const notice = renderNotice({ frontendRoot, dist, packages, frontendLicense, upstreamNotice })
  const noticePath = path.join(dist, "THIRD_PARTY_NOTICES.md")
  const licensePath = path.join(dist, "LICENSE")

  if (check) {
    if (!existsSync(noticePath) || readFileSync(noticePath, "utf8") !== notice) {
      throw new Error(`Generated notice is missing or stale: ${noticePath}`)
    }
    if (!existsSync(licensePath) || normalizeText(readFileSync(licensePath, "utf8")) !== frontendLicense) {
      throw new Error(`Frontend license is missing or stale: ${licensePath}`)
    }
    console.log(`Verified ${packages.length} production packages in ${noticePath}`)
    return
  }

  writeFileSync(noticePath, notice, "utf8")
  copyFileSync(path.join(frontendRoot, "LICENSE"), licensePath)
  console.log(`Wrote ${packages.length} production packages to ${noticePath}`)
}

if (process.argv[1] && path.resolve(process.argv[1]) === fileURLToPath(import.meta.url)) {
  main()
}
