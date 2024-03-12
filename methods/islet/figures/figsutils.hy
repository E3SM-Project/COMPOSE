(require [amb3 [*]])
(import amb3 [amb3 [*]] struct math)

(defn get-context []
  (sv c (Box)
      c.data-dir "data/mar21/"
      c.fig-dir "figs/drafts/"
      c.odes (, "rotate" "divergent" "nondivergent")
      c.cdrs (, (, "none" "none") (, "caas-node" "caas") (, "caas" "caas"))
      c.nstepfacs (, 1 5)
      c.methods (, "pcsl" "pcslu")
      c.cycs (, 1 100)
      c.timeints (, "exact" "interp")
      c.nes (, 5 10 20 40 80)
      c.nps (, 4 6 8 9 12) ;(, 4 5 6 7 8 9 10 11 12 13 16)
      c.ics (, "gau" "cos" "slo")
      c.npclrs  {4 "g" 5 "m" 6 "r" 7 "c" 8 "k" 9 "b" 10 "g" 11 "c" 12 "r" 13 "m" 16 "g"}
      c.npmarks {4 "o" 5 "x" 6 "s" 7 "x" 8 "p" 9 "+" 10 "." 11 "^" 12 "." 13 "*" 16 "."})
  c)

(defn flow-short2long [flow]
  (get {"divergent" "divergent flow"
        "nondivergent" "nondivergent flow"
        "rotate" "solid-body rotation"
        "movingvortices" "moving vortices"} flow))

(defn ic-short2long [ic]
  (get {"gau" "Gaussian hills"
        "cos" "cosine bells"
        "slo" "slotted cylinders"} ic))

(defn nes->degstrs [nes &optional [convert-all False]]
  (sv x [] xstr [])
  (for [ne nes]
    (sv deg (geton {2 15 5 6 10 3 20 "1.5" 40 "0.75" 80 "0.375" 160 "0.1875"} ne))
    (.append x ne)
    (if (none? deg)
      (do (sv deg (/ 30 ne))
          (.append xstr (.format "${:1.4f}^{{\circ}}$" deg)))
      (.append xstr (.format "${}^{{\circ}}$" deg))))
  {:ne x :deg xstr})

(defn cdr-name [short]
  (get {"caas" "CAAS-CAAS" "caas-node" "CAAS-point"} short))

;;; parse cmd, L, M, C slmmir output

(defn acc-parse [fname &optional map-nstepfac]
  (sv txt (.split (readall fname) "\n")
      bo None d {})
  (for [ln txt]
    (sv ln2 (cut ln 0 2))
    (cond [(and (= "cm" ln2) (in "cmd>" ln))
           (sv cmd ln c (parse-cmd ln :map-nstepfac map-nstepfac))]
          [(= ln2 "M ")
           (sv mp (parse-midpoint-check ln))
           (unless (none? (geton d #*(+ (cmd->key-base c)
                                        (, (:ic mp) cyc (:ne c) :M))))
             ;; handle repeated ic used for src-term OOA testing.
             (assoc mp :ic (+ (:ic mp) "-src")))
           (assoc-nested d (+ (cmd->key-base c) (, (:ic mp) cyc (:ne c) :M))
                         {:l1 (:l1 mp) :l2 (:l2 mp)})]
          [(= ln2 "C ")
           (cond [(in "cycle" ln) (sv (, - - cyc) (sscanf ln "s,s,i"))]
                 [(or (in "PASS" ln) (in "FAIL" ln))
                  (sv (, - pf) (sscanf ln "s,s"))
                  (when (and (= pf "FAIL") (!= (:mono c) "none") (<= cyc 10))
                    (prf "FAIL {}" cmd))]
                 [:else
                  (sv cl (parse-C ln))
                  (unless (none? (geton d #*(+ (cmd->key-base c)
                                               (, (:ic cl) cyc (:ne c) :C))))
                    ;; handle repeated ic used for src-term OOA testing.
                    (assoc cl :ic (+ (:ic cl) "-src")))
                  #_(print "found:" (+ (cmd->key-base c) (, (:ic cl) cyc (:ne c) :C)))
                  (assoc-nested d (+ (cmd->key-base c) (, (:ic cl) cyc (:ne c) :C))
                                cl)])]
          [(= ln2 "L ")
           (cond [(in "L file" ln)
                  (assoc-nested d (+ (cmd->key-base c)
                                     (, "cos" cyc (:ne c) :L :mixing-file))
                                (last (.split ln)))
                  ;; If nothing else has marked bo as done, do it now.
                  (unless (none? bo)
                    (assoc bo :done True))]
                 [(in "phimin" ln)
                  (sv (, - ic - l1 - l2 - li - phimin - phimax)
                      (sscanf ln "s,s,s,f,s,f,s,f,s,f,s,f"))
                  (assoc-nested d (+ (cmd->key-base c) (, ic cyc (:ne c) :Lerr))
                                {:l1 l1 :l2 l2 :li li :phimin phimin :phimax phimax})]
                 [:else
                  (sv bo (parse-bakeoff-diag bo ln (:timeint c)))])
           (when (and (not (none? bo)) (:done bo))
             (for [(, k v) (.items bo)]
               (assoc-nested d (+ (cmd->key-base c) (, "cos" cyc (:ne c) :L k)) v))
             (sv bo None))]))
  d)

;;; slmmir I/O

(defn read-slmmir-io-arrays [fname &optional beg end stride]
  (defn vis-read-array [f]
    (sv b (.read f 4))
    (when (zero? (len b)) (return None))
    (sv ndim (first (struct.unpack "@i" b))
        b (.read f (* 4 ndim))
        dims (struct.unpack (+ "@" (* "i" ndim)) b))
    (.reshape (npy.fromfile f :count (npy.prod dims)) dims))
  (svifn beg 0 end -1 stride 1)
  (sv d [])
  (with [f (open fname "rb")]
    (sv i 0 i-next beg)
    (while True
      (sv a (vis-read-array f))
      (when (none? a) (break))
      (when (= i i-next)
        (.append d a)
        (sv i-next (+ i-next stride)))
      (inc! i)))
  d)

(defn draw-slmmir-image [f &optional vmin vmax ncolor colorsym switch-halves grid]
  (svifn vmin -0.05 vmax 1.15 ncolor 24 colorsym False switch-halves True
         grid True)
  (sv (, m n) f.shape
      lon-idx (if switch-halves
                  (+ (list (range (// n 2) n)) (list (range 0 (// n 2))))
                  (s-all))
      x (* (npy.array [(/ -0.5 n) 0.25 0.5 0.75 1]) n)
      xticks [] ;(, "0" "$\pi$/2" "$\pi$" "3$\pi$/2" "$2\pi$")
      y (* (npy.array [0 0.5 1]) m)
      yticks []; (, "-$\pi$/2" "0" "$\pi$/2")
      fs 8
      colors (if colorsym
                 [(, 0 0 1) (, 1 1 1) (, 1 0 0)]
                 [(, .85 .85 .95) (, 0 0 1) (, 0 1 0) (, 1 1 0) (, 1 0 0)]))
  (if 0
      (pl.contour (get f (, (s-all) lon-idx))
                  (npy.linspace -0.05 1.15 25))
      (pl.imshow (get f (, (s-all) lon-idx))
                 (matplotlib.colors.LinearSegmentedColormap.from-list
                   "filament" colors ncolor)
                 :vmin vmin :vmax vmax))
  (pl.xlim (, (first x) (last x))) (pl.xticks x xticks :fontsize fs)
  (pl.ylim (, (first y) (last y))) (pl.yticks y yticks :fontsize fs)
  (when grid (my-grid :ls ":")))

;;; parse slmmir text output

;; (map-nstepfac nstepfac-initial ne np)
(defn parse-cmd [cmd &optional map-nstepfac]
  (sv toks (.split cmd))
  (defn int-or-none [x]
    (unless (none? x) (int x)))
  (defn get-key-val [key]
    (for [(, i t) (enumerate toks)]
      (unless (and (= (first t) "-") (= (cut t 1) key)) (continue))
      (return (get toks (inc i)))))
  (sv keys {"ode" str "ne" int "np" int "nsteps" int "prefine" int
            "mono" str "lim" str "timeint" str "method" str "pg" int-or-none}
      d {})
  (for [e (.items keys)]
    (assoc d (keyword (first e)) ((second e) (get-key-val (first e)))))
  (sv nstepfac (/ (:nsteps d) (:ne d) 6))
  (unless (none? map-nstepfac)
    (sv nstepfac (map-nstepfac nstepfac (:ne d) (:np d))))
  (assoc d :nstepfac (int nstepfac))
  (when (= (:timeint d) "exact") (assoc d :prefine 0))
  d)

(defn cmd->key-base [c]
  (sv c1 (, (:timeint c) (:ode c) (:nstepfac c) (:method c) (:mono c) (:lim c)
            (:prefine c)))
  (if (none? (:pg c))
      (+ c1 (, (:np c)))
      (+ c1 (, (:pg c) (:np c)))))

(defn parse-midpoint-check [ln]
  (sv (, - ic - l1 l2) (sscanf ln "s,s,s,f,f"))
  {:ic ic :l1 l1 :l2 l2})

(defn parse-C [ln]
  (sv (, - ic - masscons - limmin limmax - l1 l2 li - massredist massdisc)
      (sscanf ln "s,s,s,f,s,f,f,s,f,f,f,s,f,f"))
  {:ic ic :masscons masscons :limerr (, limmin limmax)
   :l1 l1 :l2 l2 :li li :massredist massredist :massdisc massdisc})

(defn parse-bakeoff-diag [d ln timeint]
  (defn parse-mixing [ln]
    (sv (, - lr - lu - lo) (sscanf (cut ln 4) "s,f,s,f,s,f"))
    {:lr lr :lu lu :lo lo})
  (defn parse-arr [ln]
    (sv toks (.split (cut ln 8)))
    (lfor t toks (float t)))
  (when (none? d) (sv d {}))
  (assoc d :done False)
  (cond [(in "   l_r" ln) (assoc d :mixing    (parse-mixing ln))]
        [(in "me l_r" ln) (assoc d :me-mixing (parse-mixing ln))]
        [(in "   thr" ln) (assoc d :thr       (parse-arr ln))]
        [(in "   fil" ln) (assoc d :fil       (parse-arr ln) :done (= timeint "exact"))]
        [(in "me fil" ln) (assoc d :me-fil    (parse-arr ln) :done True)])
  d)

;;; parse and write basis strings from search

(defn offst->Nodes [subnp offst]
  (sv nodes [])
  (for [(, i e) (enumerate subnp)]
    (sv os (get offst i))
    (.append nodes (list (range os (+ os e)))))
  nodes)

(defn Nodes->string [np nodes]
  (sv bdy 1 s (.format "{:d} {:d}" np bdy))
  (for [i (range (len nodes))]
    (+= s (.format " | {:d} {:d}:" i (len (get nodes i))))
    (for [e (get nodes i)] (+= s (.format " {:d}" e))))
  s)

(defn string->Nodes [basis-str]
  (sv toks (.split basis-str)
      np (int (first toks))
      on-bdy (int (second toks))
      nodes [] n [] ctr 0 i 3 start True)
  (assert on-bdy)
  (while (< i (len toks))
    (sv t (get toks i))
    (cond [(= t "|")
           (.append nodes n)
           (sv n [] start True)]
          [start
           (assert (= (int t) ctr))
           (inc! ctr)
           (inc! i)
           (sv start False)]
          [:else
           (.append n (int t))])
    (inc! i))
  (.append nodes n)
  (, np nodes))

(defn offset-nodal? [nodes]
  (for [n nodes]
    (sv d (npy.diff (npy.array n)))
    (when (> (npy.max d) 1) (return False)))
  True)

(defn parse-ints [s] (lfor t (.split s) (int t)))

(defn parse-search-offset-nodal-subset-line [ln]
  (sv (, - meam1 - - - - - - - wtr - npm1 npm2 npm3 - pum - - np)
      (sscanf ln "s,f,s,s,s,s,s,s,s,f,s,f,f,f,s,f,s,s,i")
      p1 (.find ln "subnp ")
      p2 (.find ln "offst ")
      subnp (parse-ints (cut ln (+ p1 5) p2))
      offst (parse-ints (cut ln (+ p2 5))))
  {:txt ln :np np :meam1 meam1 :wtr wtr :npm (, npm1 npm2 npm3) :pum pum
   :nodes (offst->Nodes subnp offst) :type :offst-nodal-subset})

(defn parse-search-nodal-subset-line [ln]
  (defn parse-nodes [s]
    (sv toks (.split s) nodess [] i 0)
    (while (< i (len toks))
      (inc! i)
      (sv nodes [])
      (.append nodess nodes)
      (while (and (< i (len toks)) (!= (nth toks i) "|"))
        (.append nodes (int (nth toks i)))
        (inc! i)))
    nodess)
  (sv (, - meam1 - - - wtr - npm1 npm2 npm3 - pum - - np)
      (sscanf ln "s,f,s,s,s,f,s,f,f,f,s,f,s,s,i")
      p1 (.find ln "subnp ")
      p2 (.find ln "nodes ")
      subnp (parse-ints (cut ln (+ p1 5) p2))
      nodes (parse-nodes (cut ln (+ p2 6))))
  {:txt ln :np np :meam1 meam1 :wtr wtr :npm (, npm1 npm2 npm3) :pum pum
   :nodes nodes :type :nodal-subset})

;; parse output from search findoffsetnodal|findnodal.
(defn parse-search-basis-line [ln]
  (if (in " offst " ln)
      (parse-search-offset-nodal-subset-line ln)
      (parse-search-nodal-subset-line ln)))

;;; JCP manuscript revision utilities

(defn jcp-context [c]
  (sv c.data-dir "data/feb24/"
      c.nps (, 4 6 8 11 13))
  c)

;; Compute ne for np given tne for np=4 such that (ne,np) has the maximal #DOF
;; <= $DOF for (tne,4).
(defn jcp-tne2ne [tne np]
  (// (* tne 3) (dec np)))

;; The floating point actual value of resolution converted to np=4.
(defn jcp-effective-tne [tne np]
  (/ (* (jcp-tne2ne tne np) (dec np)) 3))

(defn jcp-np4-ne [ne np]
  (/ (* ne (dec np)) 3))

;; Deduce nstepfac 1 or 5 from parse-cmd's internal values.
(defn jcp-nenp2nstepfac [nstepfac ne np]
  (sv nsteps (* nstepfac ne 6)
      tne (/ (* ne (dec np)) 3)
      tnstepfac (/ nsteps tne 6))
  (cond [(< tnstepfac 3) 1]
        [(< tnstepfac 8) 5]
        [:else (raisefmt "jcp-nenp2nstepfac failed: {} {} {} => {}"
                         nstepfac ne np tnstepfac)]))
