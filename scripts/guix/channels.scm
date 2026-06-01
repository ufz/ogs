(define guix-channel
  (channel
    (name 'guix)
    (url "https://codeberg.org/guix/guix.git")
    (branch "master")
    (commit "2c279760e4ea3da0848046e2495bbdf84fa7a06e")
    (introduction
     (make-channel-introduction "cdf1d7dded027019f0ebbd5d6f0147b13dfdd28d"
                                (openpgp-fingerprint
                                 "BBB0 2DDF 2CEA F6A8 0D1D  E643 A2A0 6DF2 A33A 54FA")))))

(list guix-channel)
