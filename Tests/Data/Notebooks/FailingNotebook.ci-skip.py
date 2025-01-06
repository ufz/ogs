# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.5
# ---

# %% [markdown]
# # Test Notebook assertion for testrunner
#
# A notebook test will fail with `assert False` or `raise SystemExit()`.

# %%
raise AssertionError()

# This line will not be reached but would also work:
raise SystemExit()
