(require [amb3 [*]])
(import amb3 [amb3 [*]]
        poly util
        [scipy.linalg :as linalg]
        scipy.integrate
        re math sys ctypes)

(if-main (sv amb3.*when-inp-verbosity* 1))

(defn nelem [xb yb]
  (* (dec (len xb)) (dec (len yb))))

(defn ndof [method ne np]
  (cond [(or (< method 2) (= method 5)) (* ne (** (dec np) 2))]
        [(= method 2) ne]
        [:else (raisefmt "nope")]))

(defclass Islet []
  (defn --init-- [me]
    (try (sv lib (npy.ctypeslib.load-library "libislet" ".")
             lib.power-bound.restype ctypes.c-double
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
  (defn make-space-time-op-general [me xb method np tx]
    ;; xb: eulerian mesh
    ;; tx: general target points
    ;; method: 0: natural, 1: offset nodal subset
    (sv ne (dec (len xb))
        n (* ne (dec np))
        A (npy.zeros (, (len tx) n)))
    (me.lib.make-space-time-op-for-nonuni-mesh
     (as-ctypes xb) (ctypes.c-int ne)
     (ctypes.c-int method) (ctypes.c-int np)
     (as-ctypes tx) (ctypes.c-int (len tx))
     (as-ctypes A))
    A)
  (defn make-space-time-op [me xb method np txgll]
    (sv ne (dec (len xb))
        n (* ne (dec np)))
    (assert (= (len txgll) n))
    (me.make-space-time-op-general xb method np txgll))
  (defn make-space-time-op-rons-general [me xb np subnp offst tx]
    ;; xb: eulerian mesh
    ;; tx: general target points
    ;; subnp, offst: description of offset nodal basis.
    (assert (= (len subnp) (len offst)))
    (assert (<= (len subnp) (// (inc np) 2)))
    (sv ne (dec (len xb))
        n (* ne (dec np))
        A (npy.zeros (, (len tx) n)))
    (me.lib.make-space-time-op-for-nonuni-mesh-rons
     (as-ctypes xb) (ctypes.c-int ne)
     (ctypes.c-int np) (ctypes.c-int (len subnp))
     (as-ctypes subnp) (as-ctypes offst)
     (as-ctypes tx) (ctypes.c-int (len tx))
     (as-ctypes A))
    A)
  (defn make-space-time-op-rons [me xb np subnp offst txgll]
    (sv ne (dec (len xb))
        n (* ne (dec np)))
    (assert (= (len txgll) n))
    (me.make-space-time-op-rons-general xb np subnp offst txgll))
  (defn make-space-time-op-rls-general [me xb np subnp tx]
    ;; xb: eulerian mesh
    ;; tx: general target points
    ;; subnp: description of LS basis.
    (assert (<= (len subnp) (// (inc np) 2)))
    (sv ne (dec (len xb))
        n (* ne (dec np))
        A (npy.zeros (, (len tx) n)))
    (me.lib.make-space-time-op-for-nonuni-mesh-rls
     (as-ctypes xb) (ctypes.c-int ne)
     (ctypes.c-int np) (ctypes.c-int (len subnp))
     (as-ctypes subnp)
     (as-ctypes tx) (ctypes.c-int (len tx))
     (as-ctypes A))
    A)
  (defn make-space-time-op-rls [me xb np subnp txgll]
    (sv ne (dec (len xb))
        n (* ne (dec np)))
    (assert (= (len txgll) n))
    (me.make-space-time-op-rls-general xb np subnp txgll))
  (defn make-space-time-op-nondiv2d [me xb yb method np phi-param dt]
    (assert (in method (, 0 1 2 5)))
    (sv ne (nelem xb yb)
        n (ndof method ne np)
        txgll (npy.zeros n)
        tygll (npy.zeros n)
        A (npy.zeros (, n n))
        phi-param (array-if-not phi-param))
    (me.lib.make-space-time-op-nondiv2d
     (as-ctypes xb) (ctypes.c-int (dec (len xb)))
     (as-ctypes yb) (ctypes.c-int (dec (len yb)))
     (ctypes.c-int method) (ctypes.c-int np)
     (as-ctypes phi-param) (ctypes.c-int (len phi-param))
     (ctypes.c-double dt)
     (as-ctypes txgll) (as-ctypes tygll) (as-ctypes A))
    (, txgll tygll A))
  (defn apply-space-time-op-nondiv2d [me xb yb method np phi-param dt src
                                      &optional [nstep 1] [use-caas False]]
    (assert (in method (, 0 1 2 5)))
    (sv ne (nelem xb yb)
        n (ndof method ne np)
        tgt (npy.zeros n)
        phi-param (array-if-not phi-param))
    (assert (= n (len src)))
    (me.lib.apply-space-time-op-nondiv2d
     (as-ctypes xb) (ctypes.c-int (dec (len xb)))
     (as-ctypes yb) (ctypes.c-int (dec (len yb)))
     (ctypes.c-int method) (ctypes.c-int np)
     (as-ctypes phi-param) (ctypes.c-int (len phi-param))
     (ctypes.c-double dt) (ctypes.c-int nstep) (ctypes.c-bool use-caas)
     (as-ctypes src) (as-ctypes tgt))
    tgt)
  (defn apply-space-time-density-op-div1d [me method manifold trajint np
                                           xb t0 dt nstep src]
    (assert (in method (, 0 1 2)))
    (sv ne (dec (len xb))
        n (len src)
        tgt (npy.zeros n)
        c-int ctypes.c-int c-double ctypes.c-double)
    (me.lib.apply-space-time-density-op-div1d
      (c-int method) (c-int manifold) (c-int trajint) (c-int np) (as-ctypes xb)
      (c-int ne) (c-double t0) (c-double dt) (c-int nstep) (as-ctypes src)
      (as-ctypes tgt))
    tgt)
  (defn div1d-se-interp [me xb np y xi]
    (sv yi (npy.zeros (len xi)))
    (me.lib.div1d-se-interp (as-ctypes xb) (dec (len xb)) (ctypes.c-int np)
                            (as-ctypes y) (as-ctypes xi) (ctypes.c-int (len xi))
                            (as-ctypes yi))
    yi)
  (defn div1d-cc-interp [me xb y xi]
    (sv yi (npy.zeros (len xi)))
    (me.lib.div1d-cc-interp (as-ctypes xb) (dec (len xb)) (as-ctypes y)
                            (as-ctypes xi) (ctypes.c-int (len xi)) (as-ctypes yi))
    yi)
  (defn div1d-unittest [me]
    (me.lib.div1d-unittest))
  (defn calc-power-bound [me B mu-arg]
    (when (none? me.lib) (return -1))
    (sv np (second B.shape))
    (me.lib.power-bound
     (ctypes.c-int 0)
     (as-ctypes B)
     (ctypes.c-int np)
     (ctypes.c-double mu-arg)))
  (defn apply-mixer [me xb np subnp X]
    (sv ne (dec (len xb))
        n (* ne (dec np)))
    (assert (= (first X.shape) n))
    (assert (<= subnp np))
    (me.lib.apply-mixer
     (as-ctypes xb) (ctypes.c-int ne) (ctypes.c-int np) (ctypes.c-int subnp)
     (as-ctypes X) (ctypes.c-int (second X.shape)))
    X)
  (defn get-xnodes [me method np]
    (sv xnodes (npy.zeros np))
    (me.lib.get-xnodes (ctypes.c-int method) (ctypes.c-int np) (as-ctypes xnodes))
    xnodes)
  (defn calc-xnodes-metrics-from-basis-string [me basis]
    (sv metrics (npy.zeros 3 :dtype float))
    (me.lib.calc-xnodes-metrics-from-basis-string
      (str-ctypes basis) (as-ctypes metrics))
    metrics)
  (defn calc-lebesgue-consts-from-basis-string [me basis]
    (sv metrics (npy.zeros 3 :dtype float))
    (me.lib.calc-lebesgue-consts-from-basis-string
      (str-ctypes basis) (as-ctypes metrics))
    metrics))

(defn diff [x] (- (cut x 1) (cut x 0 -1)))

(if-main
(when-inp ["islet-A" {:ne int :np int :xbopt int :txopt int :dx float}]
  (sv islet (Islet)
      xb (case/eq xbopt
                  [0
                   (npy.linspace 0 1 (inc ne))]
                  [1
                   (sort (util.np-pad 0 (npy.random.rand (dec ne)) 1))]
                  [2
                   (sv orig (npy.linspace 0 1 (inc ne))
                       perturb (npy.random.rand (inc ne))
                       perturb (* 0.75 (/ (- perturb 0.5) ne)))
                   (+ orig perturb)])
      L (- (last xb) (first xb))
      n (* ne (dec np))
      txgll (case/eq txopt
                     [0 ;; constant offset
                      (sv y (util.make-gll-mesh xb np)
                          y (+ y dx))
                      y]
                     [1 ;; random offset at element level
                      (sv y (sort (npy.random.rand (inc ne)))
                          y (* y (/ L (- (last y) (first y))))
                          y (+ dx (util.make-gll-mesh y np)))
                      y]
                     [2 ;; random u, halfway to crossing characteristics
                      (sv
                       y (util.make-gll-mesh xb np)
                       u (npy.random.rand n)
                       dt (- (* 0.5 (max (filter (fn [e] (< e 0))
                                                 (/ (diff y) (diff u))))))
                       y (+ y dx (* dt u)))
                      y]
                     [3 ;; sin(x)
                      (sv y (util.make-gll-mesh xb np :end True)
                          u (npy.sin (* 2 math.pi L y))
                          dt 0.09
                          y (+ y dx (* dt u)))
                      (cut y 0 -1)]))
  (sv delta (diff txgll)
      t-min-diff (min delta)
      t-max-diff (max delta))
  (assert (> t-min-diff 0))
  (sv delta (diff xb)
      min-diff (min delta)
      max-diff (max delta))
  (when (<= ne 3)
    (print xb)
    (print txgll)
    (print (diff txgll)))
  (for [method (, 1 0)]
    (sv A (islet.make-space-time-op xb method np txgll)
        (, lams V) (linalg.eig A :right True)
        meam1 (- (max (npy.abs lams)) 1)
        condv (util.cond-2norm V))
    (prf "{:2d} {:d} {:d} {:d} diff {:5.1f} {:4.1f} meam1: {:8.1e} condv {:8.1e} {:s}"
         np method xbopt txopt (/ t-max-diff t-min-diff) (/ max-diff min-diff)
         meam1 condv (if (> meam1 1e-14) "BAD" "g"))))

(when-inp ["test-libislet"]
  (sv islet (Islet))
  (islet.unittest)
  (sv ne 7
      np 12
      xb (util.make-cells ne :pattern 'random)
      tx (util.make-gll-mesh (util.make-cells (inc ne) :pattern 'random) np)
      A1 (islet.make-space-time-op-general xb 1 np tx)
      int (npy.dtype "int32")
      A2 (islet.make-space-time-op-rons-general
          xb np
          (npy.array [9 9 10 10 9 10] :dtype int)
          (npy.array [0 0 0 0 1 1] :dtype int) tx))
  (prf "reldif {:1.3e}" (reldif A1 A2))
  (dont with [(pl-plot (, 6 9) "fig/test-libislet")]
    (pl.subplot 2 1 1) (pl.imshow A1 :interpolation "none") (pl.colorbar)
    (pl.subplot 2 1 2) ;;(pl.imshow A2 :interpolation "none") (pl.colorbar)
    (pl.imshow (npy.log10 (npy.abs (- A2 A1))) :interpolation "none") (pl.colorbar)))

(when-inp ["test-pbc"]
  (sv islet (Islet))
  (print (islet.calc-power-bound (npy.random.randn 3 4) 0.2)))

(when-inp ["foo"]
  (sv x (npy.zeros (, 2 2))
      c (as-ctypes x))
  (print x)
  (print c)
  (sv (get x 0) 3)
  (print x)
  (print (npy.ctypeslib.as-array c))
  (print (as-ctypes (npy.array [12 9 10 10 9 10] :dtype int)))
  (print (as-ctypes (npy.array [12 9 10 10 9 10] :dtype (npy.dtype "int32")))))

(when-inp ["plot-bases" {:method int :np int :nx int :fname str}]
  ;; 0 gll_natural, 1 gll_offset_nodal_subset, 2 xnodal, 3 gll_best,
  ;; 4 uniform_offset_nodal_subset
  (sv x (npy.linspace -1 1 nx)
      islet (Islet)
      y (.transpose (islet.eval method np x)))
  (with [(pl-plot (, 4 4) fname)]
    (for [i (range np)]
      (pl.plot x (get y i) "-"))
    (my-grid)
    (axis-tight-pad :pad 0)
    (sv d 1.04)
    (pl.xlim (, (- d) d))))

(when-inp ["plot-np4" {:nx int :fname str}]
  (sv np 4
      x (npy.linspace -1 1 nx)
      islet (Islet)
      (, yn yo yb) (lfor m (, 0 1 3) (.transpose (islet.eval m np x)))
      clrs "bgrk")
  (with [(pl-plot (, 4 4) fname)]
    (for [i (range np)]
      (sv c (nth clrs i))
      (pl.plot x (get yn i) (+ c ":"))
      (pl.plot x (get yo i) (+ c "--"))
      (pl.plot x (get yb i) (+ c "-")))
    (my-grid)
    (axis-tight-pad :pad 0)
    (sv d 1.04)
    (pl.xlim (, (- d) d))))
)
