(list (channel
        (name 'guix-ogs)
        (url "https://gitlab.opengeosys.org/ogs/inf/guix-ogs.git")
        (branch "master")
        (commit "901d7522bfdb710814f26f7858fb55f57769371f"))
      (channel
        (name 'guix)
        (url "https://git.savannah.gnu.org/git/guix.git")
        (branch "master")
        (commit
          "31fe177a97bacec643180cc5bcf8805a6cb07481")
        (introduction
          (make-channel-introduction
            "cdf1d7dded027019f0ebbd5d6f0147b13dfdd28d"
            (openpgp-fingerprint
              "BBB0 2DDF 2CEA F6A8 0D1D  E643 A2A0 6DF2 A33A 54FA")))))
