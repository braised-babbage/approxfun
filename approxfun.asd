(asdf:defsystem approxfun
  :author "Erik Davis <erik@cadlag.org>"
  :license "MIT"
  :description "A toolkit for efficient computation with numerical functions via approximants."
  :serial t
  :pathname "src/"
  :components ((:file "package")
               (:file "fftw")
               (:file "chebyshev"))
  :depends-on (:cffi))
