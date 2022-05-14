(asdf:defsystem approxfun
  :author "Erik Davis <erik@cadlag.org>"
  :license "MIT"
  :description "A toolkit for efficient computation with numerical functions via approximants."
  :serial t
  :pathname "src/"
  :components ((:file "package")
               (:file "fftw")
               (:file "chebyshev")
               (:file "domain")
               (:file "approx")
               (:file "operator")
               (:file "arithmetic")
               (:file "standard-functions")
               (:file "roots")
               (:file "calculus")
               (:file "solve")
               (:file "plot"))
  :depends-on (#:alexandria
               #:cffi
	       #:interactive-gnuplot
	       #:magicl/core
	       #:magicl/ext-lapack)
  :in-order-to ((asdf:test-op (asdf:test-op #:approxfun-tests))))
