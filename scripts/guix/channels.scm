(list (channel
        (name 'guix-ogs)
        (url "https://gitlab.opengeosys.org/ogs/inf/guix-ogs.git")
        (branch "master")
        (commit "6e99eb6d84706b500e89044cd25393b4e6f3a9d6"))
      (channel
        (name 'guix)
        (url "https://git.savannah.gnu.org/git/guix.git")
        (branch "master")
        (commit
          "0c1aa354048f82023e6d2043bfecc4fd58bc9b58")
        (introduction
          (make-channel-introduction
            "cdf1d7dded027019f0ebbd5d6f0147b13dfdd28d"
            (openpgp-fingerprint
              "BBB0 2DDF 2CEA F6A8 0D1D  E643 A2A0 6DF2 A33A 54FA")))))
