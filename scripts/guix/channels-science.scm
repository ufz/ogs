(list (channel
        (name 'guix-science-nonfree)
        (url
         "https://codeberg.org/guix-science/guix-science-nonfree.git")
        (commit "65654cc6488ba787c579e3f912d507495f605092")
        (introduction
         (make-channel-introduction
          "58661b110325fd5d9b40e6f0177cc486a615817e"
          (openpgp-fingerprint
           "CA4F 8CF4 37D7 478F DA05  5FD4 4213 7701 1A37 8446"))))
      (channel
        (name 'guix-science)
        (url "https://codeberg.org/guix-science/guix-science.git")
        (commit "b4ba7cd1d7b7271b4825f033b20a4d0281796062")
        (introduction
         (make-channel-introduction
          "b1fe5aaff3ab48e798a4cce02f0212bc91f423dc"
          (openpgp-fingerprint
           "CA4F 8CF4 37D7 478F DA05  5FD4 4213 7701 1A37 8446"))))
      (channel
        (name 'guix)
        (url "https://codeberg.org/guix/guix.git")
        (branch "master")
        (commit "2c279760e4ea3da0848046e2495bbdf84fa7a06e")
        (introduction
         (make-channel-introduction
          "cdf1d7dded027019f0ebbd5d6f0147b13dfdd28d"
          (openpgp-fingerprint
           "BBB0 2DDF 2CEA F6A8 0D1D  E643 A2A0 6DF2 A33A 54FA")))))
