(require [amb3 [*]])
(import [amb3 [*]] math glob re)

(defn get-context []
  (sv c (Box)
      c.template-dir "/ccs/home/ambradl/repo/compy-goodies/slgpu/"
      c.run-dir "/ccs/home/ambradl/sl/gpu/"
      c.max-nnode 4600
      c.timers (, "main_loop" "dirk" "RK2" "tracers_compose"))
  c)

(defn get-tstep [ne] (* (/ 1024 ne) 10))
(defn get-nu [ne] (* 2.5e10 (** (/ 1024 ne) 3)))
(defn get-hv-subcyc [ne]
  (cond [(<= ne 128) 3]
        [(<= ne 512) 2]
        [:else 1]))
(defn get-remap-fac [ne] 2)
(defn get-tracer-fac [ne] 8)
(defn get-qsize [ne] 40)
(defn get-nmax [ne] (* 1 4096))
;(defn get-nmax [ne] (* 1 2048))
(defn get-nnode [ne]
  (sv nelem (* 6 (** ne 2))
      nnode (min 4096 (math.ceil (/ nelem (* 6 256)))))
  nnode)

(defn timer->sypd [nmax tstep timer]
  (/  (/ (* nmax tstep) 365)
      timer))

(defn write-files [c ne nnode-fac &optional nnode]
  ;; job: NNODE WALLTIME JOBNAME INPUTSL INPUTEUL
  ;;  nl: NE NMAX QSIZE TSTEP NUVAL HVSUB REMAPFAC TRACERFAC
  (svifn nnode -1)
  (defn make-input-name [job-name talg] (+ c.run-dir job-name "-" talg ".nl"))
  (sv nnode (if (= nnode -1) (int (* nnode-fac (get-nnode ne))) nnode)
      nmax (get-nmax ne)
      qsize (get-qsize ne)
      tstep (get-tstep ne)
      nu (get-nu ne)
      remap-fac (get-remap-fac ne)
      tracer-fac (get-tracer-fac ne)
      job-name (.format "r5-ne{}-nmax{}-qsize{}-nnode{}" ne nmax qsize nnode))
  (when (> nnode c.max-nnode) (return))
  (sed (, (, "NNODE" (str nnode)) (, "WALLTIME" "15") (, "JOBNAME" job-name)
          (, "INPUTSL" (make-input-name job-name "sl"))
          (, "INPUTEUL" (make-input-name job-name "eul")))
       (+ c.template-dir "job.sh.template")
       (+ c.run-dir (+ job-name "-job.sh")))
  (for [talg (, "eul" "sl")]
    (sv nl-template (+ "theta-" talg ".nl.template"))
    (sed (, (, "NE" (str ne)) (, "NMAX" (str nmax)) (, "QSIZE" (str qsize))
            (, "TSTEP" (str tstep)) (, "NUVAL" (.format "{:1.2e}" nu))
            (, "HVSUB" (str (get-hv-subcyc ne)))
            (, "REMAPFAC" (str remap-fac)) (, "TRACERFAC" (str tracer-fac)))
         (+ c.template-dir nl-template)
         (make-input-name job-name talg))))

(defn parse-out [c fname &optional d]
  (defn list->dict [c s]
    (assert (= (len c) (len s)))
    (sv d {})
    (for [i (range (len c))]
      (assoc d (first (nth s i))
             ((second (nth s i)) (nth c i))))
    d)
  (defn parse-fname [fname]
    (sv f (first (re.findall ".*ne(\d+)-nmax(\d+)-qsize(\d+)-nnode(\d+)" fname)))
    (list->dict f (, (, :ne int) (, :nmax int) (, :qsize int) (, :nnode int))))
  (defn parse-timer-line [ln]
    (sv (, - ngpu - ncall sum max - - - min)
        (sscanf ln "s,i,s,f,f,f,s,s,s,f"))
    {:ngpu ngpu :ncall ncall :sum sum :max max :min min})
  (svifn d {})
  (sv m (parse-fname fname)
      pat ">>>.*full")
  (for [t c.timers] (+= pat (+ "|" t ".*0\)")))
  (sv txt (grep pat fname))
  (for [ln txt]
    (cond [(= ">>>" (cut ln 0 3))
           (sv talg (second (.split ln)))]
          [:else
           (for [t c.timers]
             (when (= t (cut ln 0 (len t)))
               (when (or (and (= t "RK2") (= talg "SL"))
                         (and (= t "tracers_compose") (= talg "Eul")))
                 (raisefmt "inconsistent timers: {}" fname))
               (sv timer t
                   p (parse-timer-line (cut ln (inc (len t)))))
               (break)))
           (assoc-nested-append d (, (:ne m) (:qsize m) (:nmax m) talg
                                     timer (:nnode m))
                                p)]))
  d)

(defn parse-from-glob [c globpat]
  (sv fnames (glob.glob globpat)
      d {})
  (for [fname fnames]
    (sv d (parse-out c fname :d d)))
  d)

(defn write-table [c d &optional talgs]
  (svifn talgs (, "Eul" "SL"))
  (for [ne (sort (list (.keys d)))]
    (for [timer (, "main_loop" "advection")]
      (for [qsize (.keys (get d ne))]
        (for [nmax (.keys (get d ne qsize))]
          (sv vfirst {} first True)
          (for [talg talgs]
            (prf ">>> ne {:4d} qsize {:2d} nmax {:5d} alg {:6s} {}"
                 ne qsize nmax talg timer)
            (sv t (if (= timer "advection")
                    (if (in "SL" talg) "tracers_compose" "RK2")
                    timer)
                e (get d ne qsize nmax talg t)
                nnodes (sort (list (.keys e))))
            (for [nnode nnodes]
              (sv ps (get e nnode)
                  vs [])
              (for [p ps] (.append vs (:max p)))
              (sv v (min vs)
                  sypd (timer->sypd nmax (get-tstep ne) v)
                  speedup "")
              (if first
                (assoc vfirst nnode v)
                (sv speedup (.format "{:6.2f}" (/ (get vfirst nnode) v))))
              (prf "{:4d} {:7.2f} {:6.2f}{}" nnode v sypd speedup))
            (sv first False)))))))

;;; drivers

(when-inp ["gen"]
  (sv c (get-context))
  (for [ne (, 32 64 128 256 512 1024)
        fac (, 0.25 0.5 1 2 4)]
    (write-files c ne fac)))

(when-inp ["gen-4600"]
  (sv c (get-context))
  (for [ne (, 1024)]
    (write-files c ne 1 :nnode 4600)))

(when-inp ["table" {:globpat str}]
  (sv c (get-context)
      d (parse-from-glob c globpat))
  (write-table c d :talgs (, "Eul" "SL")))

(when-inp ["sypd" {:ne int :timer float}]
  (sv nmax (get-nmax ne)
      tstep (get-tstep ne))
  (print (timer->sypd nmax tstep timer)))

(when-inp ["fig"]
  (assoc matplotlib.rcParams "savefig.dpi" 300)
  (do (pl-require-type1-fonts))
  (sv fs 16 fsl 18)
  (defn text1 [x y dx dy data]
    (for [i (range (len x))]
      (pl.text (+ (nth x i) dx) (+ (nth y i) dy)
               (.format "{:4.2f}" (nth data i))
               :fontsize fs)))
  (defn text2 [x y dx dy data]
    (sv i (dec (len x)))
    (pl.text (+ (nth x i) dx) (+ (nth y i) dy)
             (.format "{:4.2f}" (nth data i))
             :fontsize fs))
  (defn int->str [i]
    (sv s (str i)
        s (case/eq (len s)
                   [4 (+ (first s) "," (cut s 1))]
                   [5 (+ (cut s 0 2) "," (cut s 2))]
                   [:else s]))
    s)
  (sv npa npy.array
      x (npa      [1024 2048 4096 4600])
      xfac 6
      x (* xfac x)
      sc20-y (npa [0.31 0.54 0.90 0.97])
      eul-y (npa  [0.29 0.50 0.85 0.93])
      sl-y (npa   [0.44 0.77 1.26 1.38])
      q40 (Box)
      q40.x (* xfac (npa [2048 4096 4600]))
      q40.eul-y (npa     [0.24 0.41 0.44])
      q40.sl-y (npa      [0.67 1.13 1.24])
      perf-x (* xfac (npa [1024 4600]))
      perf-y (do (sv b 0.2)
                 (npa [b (* (/ (last perf-x) (first perf-x)) b)]))
      yt (/ (npy.linspace 1 15 15) 10))
  (for [format (, "png" "pdf")]
    (with [(pl-plot (, 6 (if (= format "png") 6.6 6.4))
                    "sl-gpu-perf-032521-islet"
                    :format format)]
      (pl.plot (npy.log2 x) (npy.log2 sc20-y) "ko-" :label "SC20 data, Eulerian transport")
      (text1 (npy.log2 x) (npy.log2 sc20-y) -0.13 0.08 sc20-y)
      (pl.plot (npy.log2 x) (npy.log2 eul-y) "k.:" :label "SC20 config., Eulerian transport")
      (pl.plot (npy.log2 x) (npy.log2 sl-y) "rs-"
               :label "SC20 config., SL transport")
      (text1 (npy.log2 x) (npy.log2 sl-y) -0.13 0.08 sl-y)
      (unless (none? q40)
        (pl.plot (npy.log2 q40.x) (npy.log2 q40.eul-y) "ko--")
        (text2 (npy.log2 q40.x) (npy.log2 q40.eul-y) -0.13 0.08 q40.eul-y)
        (pl.plot (npy.log2 q40.x) (npy.log2 q40.sl-y) "rs--")
        (text2 (npy.log2 q40.x) (npy.log2 q40.sl-y) -0.03 0.03 q40.sl-y))
      (pl.plot (npy.log2 perf-x) (npy.log2 perf-y) "g:" :label "Perfect scaling")
      (pl.xticks (npy.log2 x) (lfor e x (int->str e)) :fontsize fs :rotation 45)
      (pl.yticks (npy.log2 yt) yt :fontsize fs)
      (pl.xlabel (+ "Number of Summit " (if (= xfac 6) "GPUs" "nodes"))
                 :fontsize fsl)
      (pl.ylabel "Simulated Years Per Day (SYPD)" :fontsize fsl)
      (pl.title (+ "Semi-Lagrangian tracer transport on GPU:\n"
                   "Dycore performance of SCREAM 3.25km configuration\n"
                   "solid line: 10 tracers, dashed line: 40 tracers")
                :fontsize fsl)
      (pl.legend :loc "upper left" :fontsize (dec fs) :framealpha 0)
      (my-grid)
      (sv d 0.18)
      (pl.xlim (, (- (npy.log2 (* xfac 1024)) d) (+ (npy.log2 (* xfac 4600)) d)))
      (sv d 0.27)
      (pl.ylim (, (- (npy.log (if (none? q40) 0.2 0.15)) d) (+ (npy.log 1.5) d))))))
