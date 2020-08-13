(asdf:defsystem approxfun
  :author "Erik Davis"
  :license "MIT"
  :description "A toolkit for efficient computation with numerical functions via approximants."
  :serial t
  :pathname "src/"
  :components ((:file "package")
               (:file "fftw"))
  :depends-on (:cffi))
