(list (channel
        (name 'guix-ogs)
        (url "https://gitlab.opengeosys.org/ogs/inf/guix-ogs.git")
        (branch "master")
        (commit "bbe79135fb4d6d18f3bb3ba7a3b627ffb20c5b24"))
      (channel
        (name 'guix)
        (url "https://git.savannah.gnu.org/git/guix.git")
        (branch "master")
        (commit
          "7893b32ef6ec85d51f082681a9f4fa911867afba")
        (introduction
          (make-channel-introduction
            "cdf1d7dded027019f0ebbd5d6f0147b13dfdd28d"
            (openpgp-fingerprint
              "BBB0 2DDF 2CEA F6A8 0D1D  E643 A2A0 6DF2 A33A 54FA")))))
