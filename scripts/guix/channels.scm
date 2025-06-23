(list (channel
        (name 'guix-ogs)
        (url "https://gitlab.opengeosys.org/ogs/inf/guix-ogs.git")
        (branch "master")
        (commit "8e9412330943de04e9857a75c7bf7c30904ed80e"))
      (channel
        (name 'guix-science-nonfree)
        (url "https://codeberg.org/guix-science/guix-science-nonfree.git")
        (commit "1ab2378e1d0f07ab0b591a12e741cf34a9686895")
        (introduction
         (make-channel-introduction "58661b110325fd5d9b40e6f0177cc486a615817e"
          (openpgp-fingerprint
           "CA4F 8CF4 37D7 478F DA05  5FD4 4213 7701 1A37 8446"))))
      (channel
        (name 'guix-science)
        (url "https://codeberg.org/guix-science/guix-science.git")
        (commit "b9f41a281365e0681dd78ac4c756c78a5997fa30")
        (introduction
         (make-channel-introduction "b1fe5aaff3ab48e798a4cce02f0212bc91f423dc"
          (openpgp-fingerprint
           "CA4F 8CF4 37D7 478F DA05  5FD4 4213 7701 1A37 8446"))))
      (channel
        (name 'guix)
        (url "https://codeberg.org/guix/guix.git")
        (branch "master")
        (commit "e7d73a08d569904f8a71db5b84f5fafaf0dff188")
        (introduction
         (make-channel-introduction "cdf1d7dded027019f0ebbd5d6f0147b13dfdd28d"
          (openpgp-fingerprint
           "BBB0 2DDF 2CEA F6A8 0D1D  E643 A2A0 6DF2 A33A 54FA")))))
