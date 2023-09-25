(require [amb3 [*]])
(import [amb3 [*]]
        [numpy :as npy]
        math sys)

(defn eval-lagrange-poly [x y xi]
  (setv np (len x)
        pod (pod-number? xi)
        yi (if pod 0 (npy.zeros (len xi))))
  (for [i (range np)]
    (if (= (get y i) 0) (continue))
    (setv f (if pod 1 (npy.ones (len xi))))
    (for [j (range np)]
      (if (= i j) (continue))
      (*= f (/ (- xi (get x j))
               (- (get x i) (get x j)))))
    (+= yi (* (get y i) f)))
  yi)

(defn eval-lagrange-poly-basis [x xi]
  (setv np (len x)
        pod (pod-number? xi)
        v (if pod (npy.zeros np) (npy.zeros (, np (len xi)))))
  (for [i (range np)]
    (setv f (if pod 1 (npy.ones (len xi))))
    (for [j (range np)]
      (if (= i j) (continue))
      (*= f (/ (- xi (get x j))
               (- (get x i) (get x j)))))
    (setv (get v i) f))
  v)

(defn eval-poly [coef x]
  (setv pod (pod-number? x)
        y (* (last coef) (if pod 1 (npy.ones (len x)))))
  (for [p (range (- (len coef) 2) -1 -1)]
    (setv y (+ (* x y) (get coef p))))
  y)

(defn eval-lagrange-poly-basis-derivative [x xi]
  (setv np (len x)
        pod (pod-number? xi)
        v (if pod (npy.zeros np) (npy.zeros (, np (len xi)))))
  (for [i (range np)]
    (setv f (if pod 0 (npy.zeros (len xi))))
    (for [j (range np)]
      (if (= j i) (continue))
      (sv g (if pod 1 (npy.ones (len xi))))
      (for [k (range np)]
        (if (= k i) (continue))
        (*= g (/ (if (= k j)
                   1
                   (- xi (get x k)))
                 (- (get x i) (get x k)))))
      (+= f g))
    (setv (get v i) f))
  v)

(defn eval-lagrange-poly-derivative [x y xi]
  (setv np (len x)
        pod (pod-number? xi)
        yi (if pod 0 (npy.zeros (len xi))))
  (for [i (range np)]
    (setv f (if pod 0 (npy.zeros (len xi))))
    (for [j (range np)]
      (if (= j i) (continue))
      (sv g (if pod 1 (npy.ones (len xi))))
      (for [k (range np)]
        (if (= k i) (continue))
        (*= g (/ (if (= k j)
                   1
                   (- xi (get x k)))
                 (- (get x i) (get x k)))))
      (+= f g))
    (+= yi (* (get y i) f)))
  yi)

(defn get-gll-x [np]
  (npy.array
   (case/eq np
            [1 (, 0)]
            [2 (, -1 1)]
            [3 (, -1 0 1)]
            [4 (sdo (setv oosqrt5 (/ (math.sqrt 5)))
                    (, -1 (- oosqrt5) oosqrt5 1))]
            [5 (sdo (setv sqrt3o7 (math.sqrt (/ 3 7)))
                    (, -1 (- sqrt3o7) 0 sqrt3o7 1))]
            [6 (sdo (setv e (fn [sgn]
                              (math.sqrt (+ (/ 1 3)
                                            (* sgn 2 (/ (math.sqrt 7)
                                                        21)))))
                          np6a (e 1)
                          np6b (e -1))
                    (, -1 (- np6a) (- np6b) np6b np6a 1))]
            [7 (sdo (setv e (fn [sgn]
                              (math.sqrt (/ (+ 5 (* sgn 2 (math.sqrt (/ 5 3))))
                                            11)))
                          np7a (e 1)
                          np7b (e -1))
                    (, -1 (- np7a) (- np7b) 0 np7b np7a 1))]
            [8 (sdo (setv c1 0.8717401485096066153
                          c2 0.59170018143314230214
                          c3 0.20929921790247886877)
                    (, -1 (- c1) (- c2) (- c3) c3 c2 c1 1))]
            [9 (sdo (setv c1 0.89975799541146015731
                          c2 0.67718627951073775345
                          c3 0.36311746382617815871)
                    (, -1 (- c1) (- c2) (- c3) 0 c3 c2 c1 1))]
            [10 (sdo (setv c1 0.91953390816645881383
                           c2 0.73877386510550507500
                           c3 0.47792494981044449566
                           c4 0.16527895766638702463)
                     (, -1 (- c1) (- c2) (- c3) (- c4) c4 c3 c2 c1 1))]
            [11 (sdo (setv c1 0.93400143040805913433
                           c2 0.78448347366314441862
                           c3 0.56523532699620500647
                           c4 0.29575813558693939143)
                     (, -1 (- c1) (- c2) (- c3) (- c4) 0 c4 c3 c2 c1 1))]
            [12 (sdo (setv c1 0.94489927222288222341
                           c2 0.81927932164400667835
                           c3 0.63287615303186067766
                           c4 0.39953094096534893226
                           c5 0.13655293285492755486)
                     (, -1 (- c1) (- c2) (- c3) (- c4) (- c5) c5 c4 c3 c2 c1 1))]
            [:else (raise (Exception (.format "bad np: {}" np)))])))

(defn get-gll-w [np]
  (defn reverse [coll]
    (cut coll None None -1))
  (defn expand-sym [np a]
    (if (zero? (% np 2))
      (+ a (reverse a))
      (+ a (reverse (cut a 0 -1)))))
  (when (= np 1)
    (return (, 2)))
  (npy.array
   (expand-sym
    np (case/eq
        np
        [2 (, 1)]
        [3 (, (/ 1 3) (/ 4 3))]
        [4 (, (/ 1 6) (/ 5 6))]
        [5 (, (/ 1 10) (/ 49 90) (/ 32 45))]
        [6 (sv v (math.sqrt 7))
         (, (/ 1 15) (/ (- 14 v) 30) (/ (+ 14 v) 30))]
        [7 (sv v (* 7 (math.sqrt 15)))
         (, (/ 1 21) (/ (- 124 v) 350) (/ (+ 124 v) 350) (/ 256 525))]
        [ 8 (, 0.03571428571428571429 0.21070422714350603938 0.34112269248350436476,
               0.41245879465870388157)]
        [ 9 (, 0.02777777777777777778 0.16549536156080552505 0.27453871250016173528,
               0.34642851097304634512 0.37151927437641723356)]
        [10 (, 0.02222222222222222222 0.13330599085107011113 0.22488934206312645212,
               0.29204268367968375788 0.32753976118389745666)]
        [11 (, 0.01818181818181818182 0.10961227326699486446 0.18716988178030520411,
               0.24804810426402831404 0.28687912477900808868 0.30021759545569069379)]
        [12 (, 0.01515151515151515152 0.09168451741319613067 0.15797470556437011517,
               0.21250841776102114536 0.25127560319920128029 0.27140524091069617700)]
        [:else (raisefmt "bad np: {}" np)]))))

(if-main
 (when-inp ["test-lag-basis" {:np int}]
   (for [f (, eval-lagrange-poly-basis-derivative eval-lagrange-poly-basis)]
     (setv x (npy.linspace -1 1 1000)
           v (f (get-gll-x np) x)
           v1 (f (get-gll-x np) (get x 11))))
   (expect (npy.all (= v1 (get v (, (slice None) 11)))))
   (pl.plot x (.transpose v) "-")
   (dispfig "test-lag-basis"))

 (when-inp ["plot-lagp" {:np int}]
   ;; just plot lagrange poly basis functions
   (setv x-gll (get-gll-x np)
         y-gll (npy.random.rand np)
         xi (npy.linspace -1 1 100)
         clrs "bgrcmybgrcmybgrcmy")
   (with [(pl-plot (, 6 6) "csl-plot-lagp")]
     (for [i (range np)]
       (setv y-gll (npy.zeros np)
             (get y-gll i) 1
             yi (eval-lagrange-poly x-gll y-gll xi)
             yip (eval-lagrange-poly-derivative x-gll y-gll xi)
             c (get clrs i))
       (pl.plot xi yi (+ c "-")
                xi yip (+ c "--")))
     (pl.plot x-gll (npy.zeros np) "ko")))

 (when-inp ["test-gll-w"]
   (for [np (range 1 8)]
     (sv w (get-gll-w np))
     (assert (<= (reldif 2 (sum w)) (* 1 (epsilon))))))

 (when-inp ["spacing"]
   (sv nps (list (range 4 13)) dx-min [])
   (for [np nps]
     (.append dx-min (npy.min (npy.diff (get-gll-x np)))))
   (for [i (range (len nps))]
     (sv actual (/ (first dx-min) (nth dx-min i))
         predicted (** (/ (nth nps i) (first nps)) 2))
     (prf "{:2d} {:1.3f} {:5.2f} {:4.1f} {:6.2f}" (nth nps i) (nth dx-min i)
          actual predicted (* 100 (/ predicted actual)))))
)
