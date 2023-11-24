(list (channel
        (name 'guix-ogs)
        (url "https://gitlab.opengeosys.org/ogs/inf/guix-ogs.git")
        (branch "master")
        (commit "6b51c1b7216ec1fbd6231720434c3b307bec2c0b"))
      (channel
        (name 'guix)
        (url "https://git.savannah.gnu.org/git/guix.git")
        (branch "master")
        (commit
          "0c1aa354048f82023e6d2043bfecc4fd58bc9b58")
        (introduction
          (make-channel-introduction
            "9edb3f66fd807b096b48283debdcddccfea34bad"
            (openpgp-fingerprint
              "BBB0 2DDF 2CEA F6A8 0D1D  E643 A2A0 6DF2 A33A 54FA")))))
