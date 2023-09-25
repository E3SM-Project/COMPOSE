;;; Collection of utils.

;;(require [hy.contrib.walk [let]])
(import sys copy math os re [importlib [reload]])

;;   when passing kwargs to another function like pl.plot, the dictionary should
;; be like {'option value}, not {:option value}.

(defmacro sdo [&rest code]
  "Scoped do. Just like do, but vars def'ed inside stay there."
  `((fn [] ~@code)))

(defmacro sv [&rest code]
  "The ultimate in laziness."
  `(setv ~@code))

(defmacro svb [&rest args]
  "Sub-bracketed setv, like Lisp's let."
  `(do ~@(map (fn [e]
                `(setv ~(first e) ~(last e)))
              args)))

(defn symbol [string] (HySymbol string))

(defmacro/g! svifn [&rest args]
  `(do ~@(map (fn [e]
                `(if (none? ~(first e)) (setv ~(first e) ~(last e))))
              (zip (cut args 0 (len args) 2)
                   (cut args 1 (len args) 2)))))

(defmacro/g! expect [expr &optional [answer True]]
  (setv expr-code (str expr))
  (setv answer-code (str answer))
  `(sdo (setv ~g!got ~expr)
        (setv ~g!want ~answer)
        (if (not (= ~g!got ~g!want))
          (print (.format "ERROR: {0:s} = {1:s} NOT EQUAL TO {2:s} = {3:s}"
                          ~expr-code (str ~g!got) ~answer-code (str ~g!want))))))

(defmacro/g! in-require [expr]
  (setv expr-code (str expr))
  `(sdo (setv ~g!value ~expr)
        (unless ~g!value
          (sdo (import inspect)
               (setv ~g!frame (inspect.currentframe)
                     ~g!file (. ~g!frame f-code co-filename)
                     ~g!lineno (. ~g!frame f-lineno))
               (raise (Exception (.format "IN-REQUIRE {:s} {:d}: {:s}"
                                          ~g!file ~g!lineno ~expr-code)))))))

(defmacro assert-type [v t] `(assert (= (type ~v) ~t)))

(defmacro dont [&rest code]
  "No-op."
  `((fn [])))

(defmacro raisefmt [&rest args]
  `(raise (Exception (.format ~@args))))

(defmacro interact [&rest code]
  "Block code for interactive eval, but that is silenced when the .hy file is
  run as a program."
  `(if (= --name-- "__main__")
    (dont ~@code)
    (do ~@code)))

(defmacro if-main [&rest code]
  "Block code when run with a main, but silence on import."
  `(if (= --name-- "__main__")
    (do ~@code)
    (dont ~@code)))

(defmacro prf [&rest args]
  `(print (.format ~@args)))
(defmacro prfno [&rest args]
  `(print (.format ~@args) :end ""))
(defmacro prff [fid &rest args]
  `(print (.format ~@args) :file ~fid))
(defmacro prffno [fid &rest args]
  `(print (.format ~@args) :end "" :file ~fid))

(defmacro prc [sym]
  `(print ~sym "|" (eval (read-str ~sym))))

(defmacro mapply [func &rest args]
  "Apply func. Works for keyword args, unlike apply. (Probably broken in many
  ways compared with apply, however.)"
  `(~func ~@args))

(defmacro dispatch-dict [d m]
  "Function object dispatch helper."
  `(try (get ~d ~m)
        (except [] (print "No function" (name ~m)) None)))

(defmacro/g! inc! [x &optional [increment 1]]
  `(do (setv ~x (+ ~x ~increment))
       ~x))

(defmacro/g! dec! [x &optional [increment 1]]
  `(do (setv ~x (- ~x ~increment))
       ~x))

(defmacro/g! case/test [op keyform &rest entries]
  "Case based on a test with op. Use :else as default action."
  `(do (setv ~g!keyform-result ~keyform)
       (cond ~@(map (fn [entry]
                      (+ (if (= (first entry) ':else)
                           '[True] ; If :else, use True in cond.
                           `[(~op ~g!keyform-result ~(first entry))])
                         `[~@(rest entry)]))
                    entries)
             ;; If no case matches, return None.
             [True None])))

(defmacro/g! case/eq [&rest forms]
  "Case using = for the key. Thus, the 'keylist' is not really a list, but an
  atom."
  `(case/test = ~@forms))

(defmacro/g! case/in [&rest forms]
  "Case using 'in' for the key. Thus, each keylist *must* indeed be a list."
  `(case/test in ~@forms))

(defmacro geton [&rest forms]
  "Like get, but return instead of raising KeyError on failure."
  `(try (get ~@forms)
        (except [] None)))

;; I want to switch to lfor in new code, but so old code doesn't break, provide
;; list-comp, which was removed in Hy 0.15.
(defmacro list-comp [transform range]
  `(lfor ~(first range) ~(second range) ~transform))

;; Inject the variable(s) first? or last?.
(defmacro/g! for-last [it &rest body]
  `(do
     (setv ~g!last (first (last (enumerate ~(second it)))))
     (for [(, ~g!i ~(first it)) (enumerate ~(second it))]
       (setv last? (= ~g!i ~g!last))
       ~@body)))

(defmacro/g! for-first-last [it &rest body]
  `(do
     (setv ~g!first (first (first (enumerate ~(second it))))
           ~g!last (first (last (enumerate ~(second it)))))
     (for [(, ~g!i ~(first it)) (enumerate ~(second it))]
       (setv first? (= ~g!i ~g!first))
       (setv last? (= ~g!i ~g!last))
       ~@body)))

(defclass Box []
  "A box to hold values to be written in closures."
  (defn --repr-- [me]
    (str me.--dict--)))
(defn class? [o] (= (type o) (type (Box))))
(defn has-field? [o field] (in field o.--dict--))
(defn pod-number? [n] (in (type n) (, int float)))
(defn pod-or-len [n] (if (pod-number? n) n (len n)))
(defn list? [coll] (= (type coll) list))
(defn tuple? [coll] (= (type coll) tuple))
(defn dict? [coll] (= (type coll) dict))
(defn fn? [f] (= (type f) (type (fn []))))

(defmacro/g! box-slots [&rest slots]
  "Example: (setv hi 3 bye \"so long\" b (box-slots 'hi 'bye))"
  `(do (setv ~g!box (Box))
       ~@(map (fn [s]
                ;;            handles _ vs - in names
                (setv g!field ((. (str (second s)) replace) "-" "_")
                      ;;      inject symbol directly
                      g!value (second s))
                `(assoc (. ~g!box --dict--) ~g!field ~g!value))
              slots)
       ~g!box))

(if-main
  (expect (sdo (setv hi "so long" foo-bar 3 b (box-slots 'hi 'foo-bar))
               b.foo-bar)
          3))

(defn strleq [s ref]
  (if (< (len s) (len ref))
    False
    (= (cut s 0 (len ref)) ref)))

(defn mapl [fun &rest args]
  (list (map fun #*args)))

(defn assoc-nested [d keys val]
  "Associate in a nested dict, creating new sub-dicts as needed."
  (setv dn d)
  (for [key (cut keys 0 -1)]
    (if-not (in key dn)
            (assoc dn key {}))
    (setv dn (get dn key)))
  (assoc dn (get keys -1) val))

(defn assoc-nested-append [d keys val]
  "Associate in a nested dict, creating new sub-dicts as needed. The value is
  intended to be an item to go into a list. If the list exists, append to it; if
  not, create it with the item as the only element."
  (assoc-nested d keys (+ (or (geton d #*keys) []) [val])))

;; Very nice for writing .py files from C++ with data, and then running in a
;; parse/plot program to load the data.
(defn get-python-mod [filename &optional [basedir ""]]
  "Return the module for basedir/filename."
  (defn get-module-basename [fn]
    (setv name ((. os path basename) fn))
    (get (.split name ".") 0))
  (setv name (get-module-basename filename))
  (if (not (empty? basedir))
    (+= basedir "."))
  (import importlib)
  (importlib.import-module (.format "{:s}{:s}" basedir name)))

;; Set up the path so that os.* calls know where to find .so files.
(defn set-ld-path []
  (setv (, _ stdout) (os.popen2 ". ~/sems/source.sh; env|grep LD_LIBRARY_PATH")
        paths (stdout.read)
        paths (cut paths (inc (.find paths "=")) -1)
        (get os.environ "LD_LIBRARY_PATH") paths))

;; Probably already exists, but haven't looked.
(defn unzip [coll-of-tups]
  "Unzip [(a1 b1 ...) (a2 b2 ...) ...] into [[a1 a2 ...] [b2 b2 ...] ...]."
  (if (empty? coll-of-tups)
    coll-of-tups
    (do (setv uz [])
        (for [i (range (len (get coll-of-tups 0)))]
          (.append uz (list (list-comp (get e i) [e coll-of-tups]))))
        uz)))

;; Basic cmp function that wrappers can call.
(defn cmp-fn [a b]
  (int (cond [(< a b) -1]
             [(> a b) 1]
             [:else 0])))

(defn sort! [coll]
  "Functional sort."
  (.sort coll)
  coll)

(defn sort [coll]
  "Functional sort."
  (setv c (copy.deepcopy coll))
  (sort! c))


(defn extend! [coll1 coll2]
  (.extend coll1 coll2)
  coll1)

(defn extend [coll1 coll2]
  (setv c (copy.deepcopy coll1))
  (extend! c coll2))

(defn find [e coll &optional [always-list False]]
  (setv f [])
  (for [i (range (len coll))]
    (if (= e (get coll i))
      (.append f i)))
  (if (and (not always-list) (= 1 (len f))) (first f) f))

(defn first-or-none [coll]
  (if (empty? coll) None (first coll)))

(defn safe-len [x]
  (try (len x) (except [] 1)))

(defn sscanf [str-to-parse fmt &optional split]
  """
  Kind of like sscanf. Format is like this: 'i,if,f,s', where if is for
  int(float(.)).
  """
  (setv str2type {"i" int
                  "if" (fn [x] (int (float x)))
                  "f" float
                  "s" str})
  (setv ts (list-comp (get str2type s)
                      [s (.split fmt ",")]))
  (list-comp (t s)
             [(, t s) (zip ts (.split str-to-parse split))]))

(defn split-convert [ln conversion]
  "Example:
      (split-convert \"COMPOSE> ne 24 nmax 480 qsize 5 slmpi 1\"
                     \"ssisisisi\")
   'conversion' can be shorter than (.split ln), in which case the remainder
   of line is omitted."
  (try
   (list-comp ((get {"i" int "s" identity "f" float}
                    (get conversion i)) tkn)
              [(, i tkn) (enumerate (cut (.split ln) 0 (len conversion)))])
   (except []
     (prf "split-convert failed on:\n  {:s}\nlen ln {:d} len conversion {:d}"
          ln (len (.split ln)) (len conversion)))))

(defn cc [x] (* 0.5 (+ (cut x 0 -1) (cut x 1))))

(defn mean [coll]
  (/ (sum coll) (len coll)))

(defn median [coll]
  (setv c (list (copy.deepcopy coll)))
  (.sort c)
  (if (odd? (len c))
    (get c (int (math.floor (/ (len c) 2))))
    (sdo (setv i (dec (int (/ (len c) 2))))
         (mean (cut c i (+ i 2))))))

(defn variance [coll]
  (setv mu (mean coll))
  (/ (sum (map (fn [e] (** (- e mu) 2)) coll))
     (len coll)))

(defn cumsum [coll]
  (setv c (list (copy.deepcopy coll)))
  (for [i (range 1 (len c))] (+= (get c i) (nth c (dec i))))
  c)

(defn cross-prod [x y]
  (defn one [i0 i1]
    (- (* (get x i0) (get y i1)) (* (get x i1) (get y i0))))
  [(one 1 2) (one 2 0) (one 0 1)])

(defn readall [filename]
  (with [f (open filename "r")]
    (f.read)))

(defn grep-str [pattern str]
  (import re)
  (re.findall (+ "(?:" pattern ").*") str :flags re.MULTILINE))

(defn grep [pattern filename]
  (grep-str pattern (with [f (open filename "r")] (f.read))))

(defn sed-str [pat-repls str]
  (import re)
  (for [pr pat-repls]
    (sv str (re.sub (first pr) (second pr) str :flags re.MULTILINE)))
  str)

(defn sed [pat-repls file-in file-out]
  (sv str (sed-str pat-repls (with [f (open file-in "r")] (f.read))))
  (with [f (open file-out "w")] (f.write str)))

(if-main
 (sv s (sed-str (, (,"BAR" "yes") (, "FOO" "cow"))
                "foo BAR FOO BAR\nBAR hold\nbar FOO moo"))
 (expect s "foo yes cow yes\nyes hold\nbar cow moo"))

(defn re-split-convert [converts pat ln]
  "Ex: (re-split-convert (, int float) regex ln)"
  (try
   (list-comp (convert e)
              [(, convert e) (zip converts
                                  (first (re.findall pat ln)))])
   (except [e Exception]
     (print ln))))

(defn calc-ooa [y &optional [xfac 2] [x None]]
  (defn bad [e] (or (npy.isnan e) (npy.isinf e)))
  (setv r [])
  (for [i (range (dec (len y)))]
    (when (or (bad (nth y i)) (bad (nth y (inc i))))
      (prf "calc-ooa: i {} y(i) {} y(i+1) {}" i (nth y i) (nth y (inc i))))
    (.append r (/ (- (math.log (get y i)) (math.log (get y (inc i))))
                  (if x
                    (- (math.log (get x (inc i))) (math.log (get x i)))
                    (math.log xfac)))))
  r)

(defn ooa-from-file [fname xanchor xfldno yanchor yfldno]
  "Read text from file fname, optionally scan for lines beginning with anchor,
  and read symbol fieldno, starting from 0. Return a list of OOAs."
  (sv txt (.split (readall fname) "\n")
      x [] y [])
  (for [ln txt]
    (cond [(and (>= (len ln) (len xanchor))
                (= xanchor (cut ln 0 (len xanchor))))
           (.append x (float (get (.split ln) xfldno)))]
          [(and (>= (len ln) (len yanchor))
                (= yanchor (cut ln 0 (len yanchor))))
           (.append y (float (get (.split ln) yfldno)))]))
  (, x y (ooa y :x x)))

(defn sypd [procs calls time &optional [calls-per-day 48]]
  (/ (/ calls procs calls-per-day 365)
     (/ time (* 24 3600))))

(defn single? [coll]
  (or (not (coll? coll))
      (= (len coll) 1)))

(if-main
  (expect (single? (, 'gauss)))
  (expect (single? 'gauss))
  (expect (not (single? (, 'sin 'gauss)))))

(setv *when-inp-verbosity* 2)

(defn inp [name]
  (setv b (in name sys.argv))
  (when (or (and b (> *when-inp-verbosity* 0))
            (> *when-inp-verbosity* 1)
            (= (len sys.argv) 1))
    (prf "{:s}: {:s}" (if b "DO" "av") name))
  b)

(defmacro/g! when-inp [fn-name-args &rest body]
  "Example:
     (when-inp [\"hi\" {:hi int :bye float}]
     (print hi bye))"
  (if (> (len fn-name-args) 2)
    (raise (Exception "(when-inp [fn-name &optional args] body)")))
  (unless (or (= (len fn-name-args) 1)
              (is (type (second fn-name-args)) hy.models.HyDict))
    (raise (Exception "args must be a dict")))
  (setv fn-name (first fn-name-args)
        args (if (> (len fn-name-args) 1) (second fn-name-args)))
  (defn dissect-args [args]
    (setv alist [] arg-str "")
    (for [(, i e) (enumerate (zip (cut args 0 (len args) 2)
                                  (cut args 1 (len args) 2)))]
      (.append alist
               `(setv
                 ;; grab "kw" from ":kw" and make it a symbol
                 ;;~(HySymbol (cut (first e) 2))
                 ~(HySymbol (name (first e)))
                 ;; apply type conversion
                 (try (~(second e) (get sys.argv (+ ~i 2)))
                      (except []
                        (.format "Could not parse sys.argv {:d}: {:s}"
                                 (+ ~i 2)
                                 (get sys.argv (+ ~i 2)))))))
      (+= arg-str (+ " " (name (first e)) ": " (name (second e)))))
    (, alist arg-str))
  (if args
    (do
     (setv (, alist arg-str) (dissect-args args))
     `(sdo
       (import amb3)
       (setv ~g!b (in ~fn-name sys.argv))
       (when (or (and ~g!b (> amb3.*when-inp-verbosity* 0))
                 (> amb3.*when-inp-verbosity* 1)
                 (= (len sys.argv) 1))
         (prf "{:s}: {:s}:{:s}" (if ~g!b "DO" "av")
              ~fn-name
              ~arg-str))
       (when ~g!b
         (if (< (- (len sys.argv) 2) (len ~args))
           (raise (Exception (+ "in " ~fn-name
                                " args has more entries than"
                                " are available in sys.argv"))))
         ~@alist
         ~@body)))
    `(sdo (when (inp ~fn-name) ~@body))))

(if-main
 (when-inp ["hi" {:bye int :aye float}] (print bye aye))
 (when-inp ["hello"] (print "hello")))

(defn and-coll [pred coll]
  (reduce (fn [accum e] (and accum (pred e))) coll True))

(defn or-coll [pred coll]
  (reduce (fn [accum e] (or accum (pred e))) coll False))

(defn none-in [items coll]
  (and-coll (fn [e] (not (in e coll))) items))

(defn any-in [items coll]
  (or-coll (fn [e] (in e coll)) items))

(if-main
  (expect (none-in (, "hi" "bye") "bye hello") False)
  (expect (any-in (, "hi" "bye") "bye hello"))
  (expect (none-in (, "hi" "bye") "adieu hello"))
  (expect (any-in (, "hi" "bye") "adieu hello") False)
  (expect (none-in '(1 2 3) (range 10)) False)
  (expect (none-in '(1 2 3) (range 4 10)))
  (expect (none-in (, "hi" "bye") ["bye" "hello"]) False)
  (expect (none-in (, "hi" "bye") ["adieu" "hello"])))

(defn str-ctypes [s] (ctypes.c-char-p (bytes s :encoding "ascii")))

;;; Numpy utils.

(try
 (do
  (import [numpy :as npy] ctypes)

  (defmacro getc [vs i] `(get ~vs (, (s-all) ~i)))

  (defn array-range [&rest args]
    (npy.array (list (range #*args)) :dtype int))

  (defn array-if-not [A &optional [dtype float]]
    (unless (= (type A) npy.ndarray)
      (npy.array A :dtype float)))

  (defn as-ctypes [x]
    (npy.ctypeslib.as-ctypes x))

  (defn vectorize [A] (npy.reshape A A.size))
  (defn row-vec [v] (npy.reshape v (, 1 v.size)))
  (defn col-vec [v] (npy.reshape v (, v.size 1)))
  (defn vector? [v]
    (or (and (= (len v.shape) 2) (= (min v.shape) 1))
        (= (len v.shape) 1)))
  (defn ones-row-vec [n] (npy.ones (, 1 (pod-or-len n))))
  (defn ones-col-vec [n] (npy.ones (, (pod-or-len n) 1)))

  (defn xy->XY [x y]
    (, (npy.dot (col-vec (npy.ones (len y))) (row-vec x))
       (npy.dot (col-vec y) (row-vec (npy.ones (len x))))))

  (defn sort-with-p [ai]
    "Return sorted ai and permutation array. Each entry of a must have the same
    type."
    (if (empty? ai)
      (, ai [])
      (do (setv dtype [(, (str "data") (type (get ai 0))) (, (str "p") int)]
                a (npy.array (list-comp (, e i) [(, i e) (enumerate ai)])
                            :dtype dtype))
          (setv a (npy.sort a :order (str "data")))
          (tuple (unzip a)))))

  (defn epsilon [&optional [type float]]
    (. (npy.finfo type) eps))

  (defn np-map [f a]
    (npy.array (list (map f a))))

  (defn dbg-array->np-array [a m n]
    (npy.transpose (npy.reshape (npy.array a) (, n m))))

  (defn reldif [a b &optional [norm None]]
    (if (and (pod-number? a) (pod-number? b))
      (/ (abs (- a b)) (max (abs a) (abs b)))
      (sdo (setv aa (npy.array a)
                 ba (npy.array b))
           (/ (npy.linalg.norm (- aa ba) :ord norm)
              (npy.linalg.norm aa :ord norm)))))

  (defn np-set-print []
    (setv float-format (fn [x]
                         (cond [(zero? x) (+ " ." (* " " 8))]
                               [(= x 1) (+ " 1" (* " " 8))]
                               [:else (.format "{:10.3e}" x)]))
          complex-format (fn [x]
                           (cond [(zero? x) (+ " ." (* " " 19))]
                                 [(= x 1) (+ " 1" (* " " 19))]
                                 [:else (if (zero? x.imag)
                                          (.format (+ "{:10.3e}" (* " " 11)) x.real)
                                          (.format "{:21.3e}" x))]))
          int-format (fn [x]
                       (if (zero? x)
                         (+ " ." (* " " 1))
                         (.format "{:2d} " x))))
    (npy.set-printoptions
     :precision 2
     :linewidth 1000
     :formatter {"float" float-format
                 "complexfloat" complex-format
                 "int" int-format}))

  (defn triple-read-file [fn]
    (setv triple [])
    (with [f (open fn "r")]
      (while True
        (setv ln (f.readline))
        (if (empty? ln) (break))
        (setv d (sscanf ln "i,i,f"))
        (.append triple (tuple d))))
    triple)

  (defn triple->dense [triple &optional [base 0]]
    (setv (, row col val) (unzip triple)
          A (sdo (setv d (if (= base 0) 1 0)
                       M (+ (max row) d)
                       N (+ (max col) d))
                 (npy.zeros (, M N))))
    (for [e triple]
      (setv r (get e 0)
            c (get e 1))
      (if (> base 0)
        (setv r (- r base)
              c (- c base)))
      (setv v (get e 2)
            (get A r c) v))
    A)

  (defn dense-extract-block-diag [A bs]
    (setv D (npy.zeros (npy.shape A)))
    (for [br (range (// (npy.size A 0) bs))]
      (setv r0 (* br bs)
            cs (list (range r0 (+ r0 bs))))
      (for [i (range bs)]
        (setv r (+ r0 i)
              (get (get D r) cs) (get (get A r) cs))))
    D)

  (defn dense-max-norm [A]
    (setv m1n 0)
    (for [r (range (npy.size A 0))]
      (setv r1n (as-> (get A r) it
                      (npy.abs it)
                      (npy.sum it))
            m1n (max m1n r1n)))
    m1n)

  (defn pod-number? [n]
    (in (type n) (, int float npy.int64 npy.float64)))

  (defn np-array? [a] (= (type a) npy.ndarray))

  (defn conforms? [v u]
    (and (np-array? u) (= u.shape v.shape)))

  (defn s-all [] (slice None))
  (defn s-all-rev [] (slice None None -1))

  (defn idx-arr [A rows cols]
    (get (get A rows) (, (s-all) cols)))

  (defn antidiag [v]
    (get (npy.diag (get (npy.array v) (s-all-rev))) (s-all-rev)))

  (defn get-row-2norms [vs]
    (assert (and (= (len vs.shape) 2)))
    (sv nrm (second vs.shape)
        den 0)
    (for [i (range nrm)]
      (+= den (** (getc vs i) 2)))
    (sv den (npy.sqrt den))
    den)

  (defn scale-rows! [scale vs]
    (assert (and (< (len scale.shape) 2)
                 (= (len scale) (len vs))
                 (= (len vs.shape) 2)))
    (sv dim (second vs.shape))
    (for [i (range dim)]
      (sv (getc vs i) (* scale (getc vs i))))
    vs)

  (defn normalize-rows! [vs]
    (scale-rows! (/ (get-row-2norms vs)) vs))

  )
 (except []
   (do
    (defn np-array? [a] False)
    )))

;;; Matplotlib utils.

(try
 (do
  (import matplotlib [matplotlib.pyplot :as pl])

  (defn my-grid [&optional ls]
    (svifn ls "-")
    (pl.grid True :lw 0.5 :ls ls :color (, 0.8 0.8 0.8) :zorder -1
             :which "both")
    (.set_axisbelow (pl.gca) True))

  (defn dispfig [&optional fn-prefix [format "pdf"] [tight True] [nowarn False]]
    (import warnings)
    (with [(warnings.catch-warnings)]
      (when nowarn (warnings.filterwarnings "ignore"))
      (when tight (pl.tight-layout))
      (if (or (not fn_prefix) (empty? fn-prefix))
          (pl.show)
          (pl.savefig (+ fn-prefix (+ "." format))
                      :format format :bbox-inches "tight"))))

  (defclass pl-plot []
    (defn --init-- [me figsize filename &optional [format None]
                    [tight True] [nowarn False]]
      (setv me.filename filename
            me.format (if (none? format) "pdf" format)
            me.tight tight me.nowarn nowarn)
      (pl.close)
      (pl.figure :num 1 :figsize figsize))
    (defn cleanup [me]
      (dispfig me.filename :format me.format :tight me.tight :nowarn me.nowarn))
    (defn --enter-- [me] me)
    (defn --exit-- [me &rest args] (me.cleanup))
    (defn --del-- [me]))

  ;; To get Type 1 fonts only. From
  ;;     http://nerdjusttyped.blogspot.com/2010/07/type-1-fonts-and-matplotlib-figures.html
  ;; The third one in particular really blows up rendering time, so switch this
  ;; block to False during development iterations.
  (defn pl-require-type1-fonts []
    (import matplotlib)
    (assoc matplotlib.rcParams
           "ps.useafm" True
           "pdf.use14corefonts" True
           "text.usetex" True))

  (defn imshow-matrix [A]
    (pl.imshow A :interpolation "none")
    (pl.show))

  (defn iml [A]
    (pl.imshow (npy.log10 (npy.abs A)) :interpolation "none"))

  (defn pad-lim [lim &optional [pad 0.05] [mult False]]
    (if mult
      (do (, (* (first lim) (- 1 pad))
             (* (second lim) (+ 1 pad))))
      (do (setv d (- (second lim) (first lim))
                delta (* pad d))
          (, (- (first lim) delta)
             (+ (second lim) delta)))))

  (defn axis-tight-pad [&optional [pad 0.05] [mult False]]
    (pl.axis "tight")
    (setv xl (pl.xlim) yl (pl.ylim))
    (pl.xlim (pad-lim xl pad mult))
    (pl.ylim (pad-lim yl pad mult)))

  (defn reset-colors []
    (.set-color-cycle (pl.gca) None))

  (defn good-subplot-dims [n &optional [pref-horiz False]]
    (setv d (cond [(< n 5) (case/eq n [1 (, 1 1)] [2 (, 2 1)]
                                    [3 (, 3 1)] [4 (, 2 2)])]
                  [(< n 7) (, 3 2)]
                  [(< n 10) (, 3 3)]
                  [(< n 13) (, 4 3)]
                  [(< n 17) (, 4 4)]
                  [:else (, 5 (int (math.ceil (/ n 5))))]))
    (if pref-horiz
      (, (second d) (first d))
      d))

  (defn get-linestyle-word [char]
    (get {"-" "solid" "--" "dashed" ":" "dotted" "-." "dashdot"} char))

  (defn set-tick-fontsize [fs]
    (for [ax (, "xaxis" "yaxis")]
      (sv ticks ((. (get (. (pl.gca) --dict--) ax)
                    get-major-ticks)))
      (for [tick ticks]
        ((. tick label set-fontsize) fs))))

  (defn make-reference-slope-triangle [x-span y-span slope pattern
                                       &optional [kwargs-plot None]
                                       [kwargs-text None]]
    (assert (= 2 (len x-span)))
    (assert (= 2 (len y-span)))
    (svifn kwargs-plot {})
    (svifn kwargs-text {})
    (sv (, x0 x1) x-span (, y0 y1) y-span
        dx (- x1 x0) dy (- y1 y0)
        x [x0 x1 x0 x0] y [y0 y0 y1 y0])
    (pl.plot #*[x y pattern] #**kwargs-plot)
    (pl.text #*[(+ x0 (* 0.1 dx)) (+ y0 (* 0.1 dy)) (str slope)]
             #**kwargs-text))

  (defn frame-init []
    (pl.show :block False))
  (defn frame-flip [U &optional [clim None] [pause 0.05]]
    (pl.clf)
    (pl.imshow U :interpolation "none" :origin "lower")
    (unless (none? clim)
      (pl.clim clim))
    (pl.tight-layout)
    (pl.pause pause))

  ) (except [] ))
