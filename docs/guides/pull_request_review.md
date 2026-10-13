# Pull request review

A pull request review is the process of checking whether a proposed change is correct, safe, and ready to merge.

## How to review a pull request

1. Read the PR title and description to understand the goal, scope, and any risks.
1. Review the changed files and confirm that the implementation matches the intended behavior described in the linked issue/feature specification.
1. Check that the code is clear, consistent with [project conventions](./development.md#designing-a-subcommand), and does not introduce unnecessary complexity.
1. Validate that the change contains [relevant unit tests](./development.md#unit-testing).
1. Look for edge cases, regressions, and anything that could affect reliability or maintainability.
1. Leave clear, actionable comments for the author when changes are needed.
1. Make inline suggestions where possible to enable faster implementation.
1. Approve the PR when the required/suggested changes have been correctly made and the review comments are resolved; request changes if important issues remain.
1. Inform the reviewee if you will not be available for review within a week and give an estimate to when your review could be expected. If you will not have time for a review, or if you deem that you do not possess the right expertise, decline the review.

## Good review habits

- Keep feedback specific and constructive.
- Focus on correctness, readability, and maintainability.
- Avoid broad scope changes unrelated to the PR - if desired they can be added as an issue for later implementation.
- Confirm that tests or validation evidence are included before merging.

This keeps reviews focused, consistent, and easier for contributors to act on.