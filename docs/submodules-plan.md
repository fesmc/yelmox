# Submodule conversion plan for yelmox

Future-reference plan for converting the four sibling-directory symlinks in the yelmox root (`fesm-utils`, `FastIsostasy`, `rembo1`, `yelmo`) into proper git submodules. Not yet executed.

## Goal

Have yelmox pin specific commits of each component repo and clone them automatically as submodules under the existing paths, replacing the current symlinks to sibling clones under `/Users/alrobi001/models/`.

## Current state

yelmox root contains symlinks to sibling clones:

| symlink        | remote                                               |
|----------------|------------------------------------------------------|
| `fesm-utils`   | `git@github.com:fesmc/fesm-utils.git`                |
| `FastIsostasy` | `git@github.com:palma-ice/FastIsostasy.git`          |
| `rembo1`       | `git@github.com:alex-robinson/rembo1.git`            |
| `yelmo`        | `git@github.com:fesmc/yelmo.git`                     |
| `ice_data`     | (leave as symlink — not in scope)                    |

These names are listed in `.gitignore`. The Makefile already uses inner paths like `fesm-utils/utils/include-serial`, `FastIsostasy/libisostasy/include`, `yelmo/libyelmo/include`, `rembo1/librembo/include`, which will be preserved by submodules.

## Decisions made

- **URLs:** SSH form, matching current remotes. For `yelmo` use `git@github.com:fesmc/yelmo.git` (not the `palma-ice` mirror).
- **Branch to track:** `main` for all four (default branch on each origin).
- **`.gitmodules` form:** pin to commit only — *no* `branch = …` line. Bumping the pin is an explicit `git add <submodule>` in yelmox.
- **Initial pin:** tip of `origin/main` at the time of conversion. Re-fetch and re-check tips before running.
- **`ice_data`:** untouched.

## Steps to execute

From the yelmox root:

```bash
# 1. Remove symlinks (no data lost — sibling clones at /Users/alrobi001/models/<name> remain)
rm fesm-utils FastIsostasy rembo1 yelmo

# 2. Drop the four entries from .gitignore (paths will now be tracked as submodules)
#    Keep: ice_data, libyelmox, unit_tests, libs/lis*, isostasy_data, etc.

# 3. Add each submodule (clones fresh into yelmox/<path>)
git submodule add git@github.com:fesmc/fesm-utils.git       fesm-utils
git submodule add git@github.com:palma-ice/FastIsostasy.git FastIsostasy
git submodule add git@github.com:alex-robinson/rembo1.git   rembo1
git submodule add git@github.com:fesmc/yelmo.git            yelmo

# 4. Verify .gitmodules + each submodule pin, then single commit:
git commit -m "add submodules: fesm-utils, FastIsostasy, rembo1, yelmo"
```

## Post-conversion: build steps still required

Each submodule contributes only its source. The build artifacts the Makefile depends on (`fesm-utils/fftw-serial/`, `fesm-utils/lis-serial/`, `FastIsostasy/libisostasy/include/`, `yelmo/libyelmo/include/`, `rembo1/librembo/include/`) are gitignored inside each submodule and must be built per the component's own install/build process — the same as today.

## Recommended config (one-time, global)

```bash
git config --global push.recurseSubmodules check
```

So `git push` in yelmox refuses if any submodule has unpushed commits at its pinned sha. Catches the easy mistake of pushing yelmox with a yelmo pin nobody else can fetch.

## Development workflow on a submodule

```bash
cd yelmox/yelmo
git checkout main                   # ALWAYS — fresh clones land in detached HEAD
# edit, build via yelmox Makefile, iterate
git commit -am "..."
git push                            # to fesmc/yelmo

cd ..                               # back to yelmox root
git status                          # shows: modified: yelmo (new commits)
git add yelmo                       # bump the pin
git commit -m "bump yelmo to <sha>: <reason>"
```

After a fresh `git clone --recurse-submodules` of yelmox elsewhere, put all submodules on their branches in one shot:

```bash
git submodule foreach 'git checkout main || git checkout master'
```

### Coordinated cross-submodule changes

If a yelmox change needs a matching yelmo change: push the yelmo branch first, then in yelmox commit the gitlink bump *together with* the yelmox-side code that uses the new API — that single yelmox commit becomes the atomic unit. `push.recurseSubmodules = check` enforces the ordering.

### Relationship to existing sibling clones

After conversion there are **two clones** of (e.g.) yelmo on disk: the original at `/Users/alrobi001/models/yelmo` and the submodule at `/Users/alrobi001/models/yelmox/yelmo`. They drift unless synced via the remote. Three options at conversion time:

- **(a)** Abandon the siblings; develop only inside the submodule. Simplest mental model.
- **(b)** Keep the siblings as the primary workspace; push from sibling and `git pull` inside the submodule. Two extra commands per sync but reuses one dev tree across multiple yelmox checkouts/worktrees.
- **(c)** Point the submodule at the sibling locally with `git config submodule.<name>.url /Users/alrobi001/models/<name>` (then `git submodule sync`). Submodule fetches from the local path; `.gitmodules` still points at GitHub for everyone else. Useful for one-yelmo-feeds-many-yelmox setups.

## Open follow-ups (deferred)

- `fesm-utils` local clone has uncommitted edits on `alex-dev` to `install.sh`, `utils/src/esm.f90`, `utils/test/test_esm.nml`, plus untracked `fesm-utils/`, `utils/src/esm_base.f90`, `utils/src/esm_extra.f90`. These are not on the origin/main tip and will not be in the initial pin. Decide before conversion whether to land them on `main` first.
- After conversion, the `fesm-utils` + `FastIsostasy` `alex-dev` branches still exist and are ahead of/behind `main` in various ways — decide separately whether to retire those branches or keep them as long-lived dev lines.

## Pins captured at planning time (2026-05-21)

For reference only; re-fetch before executing.

| submodule    | sha       | message                                                |
|--------------|-----------|--------------------------------------------------------|
| fesm-utils   | f7a3917   | Merge PR #4 from fesmc/alex-dev                        |
| FastIsostasy | 2969564   | Merge PR #46 from palma-ice/adaptive-dt                |
| rembo1       | aa999a9   | Updated runme.                                         |
| yelmo        | b3e233f   | docs: move Margin-front mass balance ...               |
