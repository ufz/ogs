(append (list (channel
                (name 'guix-science-nonfree)
                (url
                 "https://codeberg.org/guix-science/guix-science-nonfree.git")
                (commit "7bf39ac5fbe36905b87a0632945a29c7c67e31de")
                (introduction
                 (make-channel-introduction
                  "58661b110325fd5d9b40e6f0177cc486a615817e"
                  (openpgp-fingerprint
                   "CA4F 8CF4 37D7 478F DA05  5FD4 4213 7701 1A37 8446"))))
              (channel
                (name 'guix-science)
                (url "https://codeberg.org/guix-science/guix-science.git")
                (commit "4854c8992cdfeee4e537c8f36af04c110116a5c4")
                (introduction
                 (make-channel-introduction
                  "b1fe5aaff3ab48e798a4cce02f0212bc91f423dc"
                  (openpgp-fingerprint
                   "CA4F 8CF4 37D7 478F DA05  5FD4 4213 7701 1A37 8446")))))
        (load "channels.scm"))
