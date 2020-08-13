(fiasco:define-test-package #:approxfun-tests
    (:use #:cl #:approxfun)
  (:export #:run-approxfun-tests))

(in-package #:approxfun-tests)

(defun run-approxfun-tests ()
  (run-package-tests :package ':approxfun-tests
                     :verbose t
                     :describe-failures t
                     :interactive nil))
