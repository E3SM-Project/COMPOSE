(require [amb3 [*]])
(import amb3 [amb3 [*]]
        [scipy.linalg :as linalg]
        scipy.integrate
        re math sys ctypes)

(defn nelem [xb yb]
  (* (dec (len xb)) (dec (len yb))))

(defn ndof [method ne np]
  (cond [(< method 2) (* ne (** (dec np) 2))]
        [(= method 2) ne]
        [:else (raisefmt "nope")]))

(defclass Islet []
  (defn --init-- [me]
    (try (sv lib (npy.ctypeslib.load-library "libislet" ".")
             me.lib lib)
         (except [e [Exception]]
           (print e)
           (sv me.lib None))))
  (defn unittest [me]
    (me.lib.islet-unittest))
  (defn eval [me method np x]
    (sv c-int ctypes.c-int
        y (npy.zeros (, (len x) np)))
    (me.lib.eval-interpolant (c-int method) (c-int np) (c-int (len x))
                             (as-ctypes x) (as-ctypes y))
    y)
  (defn get-xnodes [me method np]
    (sv xnodes (npy.zeros np))
    (me.lib.get-xnodes (ctypes.c-int method) (ctypes.c-int np) (as-ctypes xnodes))
    xnodes)
  (defn calc-xnodes-metrics-from-basis-string [me basis]
    (sv metrics (npy.zeros 3 :dtype float))
    (me.lib.calc-xnodes-metrics-from-basis-string
      (str-ctypes basis) (as-ctypes metrics))
    metrics))

(defn diff [x] (- (cut x 1) (cut x 0 -1)))

