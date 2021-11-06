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
               (:file "plot"))
  :depends-on (#:alexandria #:cffi #:magicl)
  :in-order-to ((asdf:test-op (asdf:test-op #:approxfun-tests))))
