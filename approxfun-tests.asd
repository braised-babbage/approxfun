(asdf:defsystem approxfun-tests
  :author "Erik Davis <erik@cadlag.org>"
  :license "MIT"
  :description "Tests for approxfun."
  :perform (asdf:test-op (o s)
                         (uiop:symbol-call 'approxfun-tests
                                           '#:run-approxfun-tests))
  :serial t
  :pathname "tests/"
  :components ((:file "test"))
  :depends-on (#:approxfun #:fiasco))
