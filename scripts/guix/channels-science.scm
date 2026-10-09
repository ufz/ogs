(load (string-append (dirname (current-filename)) "/channels.scm"))

(list (channel
        (name 'guix-science-nonfree)
        (url
         "https://codeberg.org/guix-science/guix-science-nonfree.git")
        (commit "b9cc745a1d756cf94de355a12d5742c71230fe7d")
        (introduction
         (make-channel-introduction
          "58661b110325fd5d9b40e6f0177cc486a615817e"
          (openpgp-fingerprint
           "CA4F 8CF4 37D7 478F DA05  5FD4 4213 7701 1A37 8446"))))
      (channel
        (name 'guix-science)
        (url "https://codeberg.org/guix-science/guix-science.git")
        (commit "26bc52b6d3e502e84b2c4b8c1033f3e1b5d5bfa8")
        (introduction
         (make-channel-introduction
          "b1fe5aaff3ab48e798a4cce02f0212bc91f423dc"
          (openpgp-fingerprint
           "CA4F 8CF4 37D7 478F DA05  5FD4 4213 7701 1A37 8446"))))
      guix-channel)
