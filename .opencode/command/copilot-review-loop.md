---
description: Run the Copilot review loop on a PR until it converges
agent: build
---

# Copilot review loop

You are running the Copilot review loop on PR **$ARGUMENTS**. The loop is the contract; one round is not enough. Resolve `<owner>`, `<repo>`, and `<pr>` from `$ARGUMENTS` (accept a number, a PR URL, or `owner/repo#N`; default to the current repo via `gh repo view --json nameWithOwner` if only a number is given).

## Requesting a review

The only request form that works is the literal-bot-suffix form of the REST API:

```bash
gh api --method POST \
  /repos/<owner>/<repo>/pulls/<pr>/requested_reviewers \
  -f 'reviewers[]=copilot-pull-request-reviewer[bot]'
```

`copilot-pull-request-reviewer` (without the `[bot]` suffix) returns `422 Reviews may only be requested from collaborators`. `Copilot` returns 200 but silently drops the request. The `requestReviews` GraphQL mutation returns `NOT_FOUND` for the bot's node id. None of those work. The literal `[bot]` suffix on the reviewer name is the only path; the successful response shows `requested_reviewers: [{login: "Copilot", type: "Bot"}]` and the bot posts a review within a few minutes.

## The loop

1. **Request a review** using the command above.
2. **Watch for the review.** Poll BOTH endpoints in tandem: `gh api repos/<owner>/<repo>/pulls/<pr>/reviews` for the review count, and `gh api repos/<owner>/<repo>/pulls/<pr>/requested_reviewers` for the bot's state. After the request, `requested_reviewers` shows `[{login: "Copilot", type: "Bot"}]` (the bot has the request queued); it clears back to `[]` when the bot finishes processing. Two outcomes from there:
   - **A new review IS posted** (`reviews` count rises): continue to step 3.
   - **`requested_reviewers` cleared and no new review was posted**: the bot processed the request and had nothing to add. Stop the loop now. **Do not run a confirmation pass**; the cleared `requested_reviewers` is itself the terminating signal.

   Five minutes is a reasonable cap for the bot to process a request. If `requested_reviewers` has not cleared after that, request again.
3. **Read every comment in the new review.** For each thread, either fix or reject:
   - **Fix** if the finding is correct. Keep commits focused — one concern per commit; restructures, renames, and reorganisations split into their own commits.
   - **Reject** if the finding is a false positive, a misreading of the code, a misinterpretation of a deliberate decision (e.g. the documented SI-unit contract, an intentional `NaN`/`Inf` domain edge, a cited formula), or otherwise does not warrant a change. Do NOT silently apply suggestions that would degrade the code, contradict documented design, or revert deliberate work. A rejected finding gets a reply that names the reason in one sentence.
4. **Reply to each comment** with `gh api repos/<owner>/<repo>/pulls/<pr>/comments -f body=... -F in_reply_to=<comment_database_id>` and **resolve each thread** via the GraphQL `resolveReviewThread` mutation against the thread node id.
5. **Push the fix commits** before requesting the next round so Copilot reviews the new state, not the stale one.
6. **Request another review** and repeat from step 2.
7. **Stop the loop** when one of these holds: Copilot processed a request without posting a new review (`requested_reviewers` cycled `[Copilot] → []` and the reviews count did not rise); or the latest posted review surfaces no new actionable findings (review body says "no issues", or every comment is a duplicate of one already resolved, or every comment is a false positive already rejected with reasoning). The loop terminates on Copilot's behaviour, not on a fixed count, and a single round meeting either condition is enough; no confirmation pass.

Begin now with step 1 for PR $ARGUMENTS.
