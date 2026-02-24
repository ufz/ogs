(define guix-channel
  (channel
    (name 'guix)
    (url "https://codeberg.org/guix/guix.git")
    (branch "master")
    (commit "4b5932da6d0bf3f8987643f8c1271d1816cba698")
    (introduction
     (make-channel-introduction "cdf1d7dded027019f0ebbd5d6f0147b13dfdd28d"
                                (openpgp-fingerprint
                                 "BBB0 2DDF 2CEA F6A8 0D1D  E643 A2A0 6DF2 A33A 54FA")))))

(list guix-channel)
