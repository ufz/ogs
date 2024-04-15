(list (channel
        (name 'guix-ogs)
        (url "https://gitlab.opengeosys.org/ogs/inf/guix-ogs.git")
        (branch "master")
        (commit "151fe23f8610b85a84dbaf6888a5afc1256cbb10"))
      (channel
        (name 'guix)
        (url "https://git.savannah.gnu.org/git/guix.git")
        (branch "master")
        (commit
          "6cb181c07f83dfdeae1882208941086f3717a165")
        (introduction
          (make-channel-introduction
            "cdf1d7dded027019f0ebbd5d6f0147b13dfdd28d"
            (openpgp-fingerprint
              "BBB0 2DDF 2CEA F6A8 0D1D  E643 A2A0 6DF2 A33A 54FA")))))
