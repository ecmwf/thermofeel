---
description: Address all unresolved review comments on a GitHub PR and push fixes
agent: build
---

# Address Pull Request Comments

Review and resolve all feedback on a GitHub pull request.

## Arguments

Provide the PR number as argument, e.g. `/address-pr-comments 32`
If no number given, detect the current branch's open PR.

## Process

### 1. Fetch Reviews
First, fetch all comments in **unresolved** threads and reviews:
```bash
gh api graphql -f query='
{
  repository(owner: "ecmwf", name: "thermofeel") {
    pullRequest(number: PR_NUMBER) {
      reviewThreads(first: 50) {
        nodes {
          isResolved
          comments(first: 1) {
            nodes { path line originalLine body }
          }
        }
      }
    }
  }
}' | jq -r '.data.repository.pullRequest.reviewThreads.nodes[] | select(.isResolved == false) | "File: \(.comments.nodes[0].path)\nLine: \(.comments.nodes[0].line // .comments.nodes[0].originalLine)\nComment: \(.comments.nodes[0].body)\n---"'
```
This provides the actual inline review comments which `gh pr view` doesn't show properly, and filters out resolved threads.

Also obtain reviews and root-level comments with
```
gh api repos/{owner}/{repo}/pulls/<NUMBER>/comments
gh api repos/{owner}/{repo}/pulls/<NUMBER>/reviews
```

### 2. Analyze Each Comment

For each review comment or inline comment:
- Read the full context of the code being discussed
- Understand the reviewer's concern in light of this project's philosophy:
  - **Scientific correctness first** — every formula matches its cited
    reference; validate against published worked examples, not a self-generated
    number
  - **The SI-unit in/out contract** (Kelvin, m/s, %, hPa) is upheld at the
    public boundary (see `plans/DESIGN.md`)
  - **This is a published library** — public signatures in
    `thermofeel/thermofeel.py` and the `helpers.py` converters are a contract;
    behaviour changes follow SemVer and are recorded in `CHANGELOG.md`
  - Proper fixes over workarounds; no suppressed lints/tests
- If the intent is ambiguous, ask the user for clarification before acting

### 3. Address Issues

For each actionable comment:
- Fix the code as requested (or propose a better alternative with justification)
- Update documentation (`docs/`, docstrings) if the change affects public API or
  behaviour, and `CHANGELOG.md` for any user-facing change
- Add or update tests if the comment identified a gap (scalar and/or array
  suite, plus `tests/test_robustness.py` for edge/`NaN` behaviour)

### 4. Push and Verify

- Commit changes with a clear message referencing the review feedback
- Push to update the PR
- Wait for CI to finish:
  ```bash
  gh pr checks <NUMBER> --watch
  ```
- If CI fails, fix and push again

### 5. Iterate

Continue the loop until:
- All review comments have been addressed
- CI is green
- Summarize what was changed in response to each comment

### 6. Report

Provide a summary mapping each reviewer comment to the action taken.
