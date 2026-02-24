(append (list (channel
                (name 'guix-science-nonfree)
                (url
                 "https://codeberg.org/guix-science/guix-science-nonfree.git")
                (commit "4afed2444866b85b44ba9b0dceb874917ca07df7")
                (introduction
                 (make-channel-introduction
                  "58661b110325fd5d9b40e6f0177cc486a615817e"
                  (openpgp-fingerprint
                   "CA4F 8CF4 37D7 478F DA05  5FD4 4213 7701 1A37 8446"))))
              (channel
                (name 'guix-science)
                (url "https://codeberg.org/guix-science/guix-science.git")
                (commit "a91b0e7c031e7da82d514bcba4b46fde070d78a5")
                (introduction
                 (make-channel-introduction
                  "b1fe5aaff3ab48e798a4cce02f0212bc91f423dc"
                  (openpgp-fingerprint
                   "CA4F 8CF4 37D7 478F DA05  5FD4 4213 7701 1A37 8446")))))
        (load "channels.scm"))
