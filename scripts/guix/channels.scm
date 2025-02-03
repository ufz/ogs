(list (channel
        (name 'guix-ogs)
        (url "https://gitlab.opengeosys.org/ogs/inf/guix-ogs.git")
        (branch "master")
        (commit "559998fedc0cadd67f5a6ab4138e928978eb9302"))
      (channel
        (name 'guix-science-nonfree)
        (url "https://codeberg.org/guix-science/guix-science-nonfree.git")
        (commit "7bb0689a9b02b96b75a527f43a6806872f1227cb")
        (introduction
         (make-channel-introduction "58661b110325fd5d9b40e6f0177cc486a615817e"
          (openpgp-fingerprint
           "CA4F 8CF4 37D7 478F DA05  5FD4 4213 7701 1A37 8446"))))
      (channel
        (name 'guix-science)
        (url "https://codeberg.org/guix-science/guix-science.git")
        (commit "85e684db279109050e2b20e6ae62b2b5d3bf5eec")
        (introduction
         (make-channel-introduction "b1fe5aaff3ab48e798a4cce02f0212bc91f423dc"
          (openpgp-fingerprint
           "CA4F 8CF4 37D7 478F DA05  5FD4 4213 7701 1A37 8446"))))
      (channel
        (name 'guix)
        (url "https://git.savannah.gnu.org/git/guix.git")
        (branch "master")
        (commit "29d133063010cc4f7c40aadb3124364eaf276797")
        (introduction
         (make-channel-introduction "cdf1d7dded027019f0ebbd5d6f0147b13dfdd28d"
          (openpgp-fingerprint
           "BBB0 2DDF 2CEA F6A8 0D1D  E643 A2A0 6DF2 A33A 54FA")))))
