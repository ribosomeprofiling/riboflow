// Groovy helpers shared by the workflow scripts and process script blocks.
//
// NOTE: conf/*.config closures CANNOT see this class (a bare `Utils` there resolves
// to a groovy.util.ConfigObject). Anything a config closure needs lives as a
// file-local closure at the top of conf/modules.config; everything param-derived
// that a *module* needs is resolved here from `params` + an `ext.route` string, so
// the config never has to spell the defaults out.

class Utils {

    // ── booleans ────────────────────────────────────────────────────────────
    // Params arrive as Boolean from YAML/config but as String from the CLI and
    // from JSON params files, so `x == true` / `x != false` mis-evaluates "true"
    // and "false". Always resolve through here.
    static boolean as_bool(Object v, boolean dflt) {
        if (v == null) return dflt
        if (v instanceof Boolean) return v
        def s = v.toString().trim().toLowerCase()
        if (s in ['true', 'yes', 'on', '1'])  return true
        if (s in ['false', 'no', 'off', '0', '']) return false
        throw new IllegalArgumentException("expected a boolean but got `${v}`")
    }

    // ── dedup method ────────────────────────────────────────────────────────
    // Normalise the dedup_method param, honouring the legacy boolean `deduplicate`
    // flag. (RiboFlow.groovy:39-62)
    static String dedup_method(String dedup_arg, String dedup_old) {
        def valid_methods = ['position', 'umicollapse', 'none']
        def dedup_param = dedup_arg.toLowerCase()
        if (dedup_param != 'none') {
            if (dedup_param in valid_methods) return dedup_param
            System.err.println('Invalid deduplication method ' + dedup_param +
                               ' . Valid methods are: ' + valid_methods.join(','))
            System.exit(1)
        }
        return as_bool(dedup_old, false) ? 'position' : 'none'
    }

    static String resolve_dedup_method(Map params) {
        return dedup_method((params.get('dedup_method', 'none') ?: 'none').toString(),
                            (params.get('deduplicate', false)).toString())
    }

    static String resolve_rnaseq_dedup_method(Map params) {
        return dedup_method((params.rnaseq?.get('dedup_method', 'none') ?: 'none').toString(), 'false')
    }

    static boolean do_tx_dedup(Map params) {
        def m = resolve_dedup_method(params)
        return (m == 'umicollapse' || m == 'position') &&
               as_bool((params.star ?: [:]).output_transcriptome_bam, false)
    }

    // ── resources ───────────────────────────────────────────────────────────
    // Per-thread memory budget for `samtools sort`. (RiboFlow.groovy:64-79)
    static int samtools_sort_mem_per_thread_mb(task) {
        int sort_threads = Math.min(task.cpus as int, 8)
        int est = (int) (task.memory.toMega() * 0.7 / sort_threads)
        return Math.min(768, Math.max(64, est))
    }

    // Memory budget (MB) for an external `sort`: 60 % of the task allocation.
    static int sort_mem_mb(task) {
        return task.memory ? Math.max(256, (int) (task.memory.toMega() * 0.6)) : 1024
    }

    // ════════════════════════════════════════════════════════════════════════
    //  Alignment routes and their post-alignment quality filter
    //
    //  Four routes, each with its own params block:
    //    genome                → params.genome
    //    transcriptome         → params.transcriptome
    //    rnaseq_genome         → params.rnaseq.genome (flat rnaseq.* keys honoured)
    //    rnaseq_transcriptome  → params.rnaseq.transcriptome
    //  plus 'star_transcriptome' (the STAR --quantMode TranscriptomeSAM projection),
    //  which shares the transcriptome keys.
    //
    //  Per-route keys:
    //    unique_only               genome routes only. true ⇒ keep uniquely-mapping
    //                              reads (samtools view -q 255; STAR gives unique
    //                              reads MAPQ 255) and emit the unique-only stats
    //                              layout. false ⇒ apply mapping_quality_cutoff and
    //                              emit the unique/multi breakdown. When absent the
    //                              mode is inferred from the cutoff (255 ⇒ unique on
    //                              ribo-seq; >0 ⇒ unique on RNA-seq) with a warning.
    //    mapping_quality_cutoff    samtools view -q (ignored while unique_only: true)
    //    samtools_filter_arguments any other `samtools view` record-selection args
    //    samtools_count_arguments  the flags behind the stats counters, see below
    // ════════════════════════════════════════════════════════════════════════

    static final List<String> ROUTES        = ['genome', 'transcriptome', 'rnaseq_genome', 'rnaseq_transcriptome']
    static final List<String> GENOME_ROUTES = ['genome', 'rnaseq_genome']
    static final int STAR_UNIQUE_MAPQ       = 255

    static boolean is_genome_route(String route) { return route in GENOME_ROUTES }

    static Map route_params(Map params, String route) {
        switch (route) {
            case 'genome':
                return ((params.genome ?: [:]) as Map)
            case 'transcriptome':
            case 'star_transcriptome':
                return ((params.transcriptome ?: [:]) as Map)
            case 'rnaseq_genome':
                def rg = ((params.rnaseq ?: [:]) as Map)
                def g  = ((rg.genome ?: [:]) as Map)
                // Flat legacy shape: rnaseq.mapping_quality_cutoff / rnaseq.samtools_filter_arguments
                def merged = [:]
                ['mapping_quality_cutoff', 'samtools_filter_arguments', 'unique_only', 'samtools_count_arguments'].each { k ->
                    if (rg.containsKey(k) && !g.containsKey(k)) merged[k] = rg[k]
                }
                return merged + g
            case 'rnaseq_transcriptome':
                return ((((params.rnaseq ?: [:]) as Map).transcriptome ?: [:]) as Map)
            default:
                throw new IllegalArgumentException("unknown alignment route `${route}`")
        }
    }

    static int default_mapq(String route) {
        switch (route) {
            case 'genome':        return STAR_UNIQUE_MAPQ
            case 'rnaseq_genome': return 4
            default:              return 10
        }
    }

    // `!= null`, NEVER Elvis: a deliberate cutoff of 0 (keep multimappers) must not
    // fall through to the default.
    static int route_mapq_cutoff(Map params, String route) {
        def v = route_params(params, route).mapping_quality_cutoff
        return (v != null) ? (v as int) : default_mapq(route)
    }

    static boolean route_unique_only_declared(Map params, String route) {
        return is_genome_route(route) && route_params(params, route).unique_only != null
    }

    static boolean route_unique_only(Map params, String route) {
        if (!is_genome_route(route)) return false
        def v = route_params(params, route).unique_only
        if (v != null) return as_bool(v, true)
        def mapq = route_mapq_cutoff(params, route)
        return (route == 'genome') ? (mapq >= STAR_UNIQUE_MAPQ) : (mapq > 0)
    }

    // The -q that SAMTOOLS_QPASS actually applies.
    static int route_effective_mapq(Map params, String route) {
        return route_unique_only(params, route) ? STAR_UNIQUE_MAPQ : route_mapq_cutoff(params, route)
    }

    // Default record filter: drop unmapped(4) + secondary(256) + supplementary(2048).
    static final String DEFAULT_SAMTOOLS_FILTER_ARGS = '-F 2308'

    // `!= null`: `samtools_filter_arguments: ""` (filter nothing) is honoured as written.
    static String route_filter_args(Map params, String route) {
        def v = route_params(params, route).samtools_filter_arguments
        return (v != null) ? v.toString() : DEFAULT_SAMTOOLS_FILTER_ARGS
    }

    // ── samtools_count_arguments ────────────────────────────────────────────
    // Flags behind the four `samtools view -c` stats counters. Paired-end keys exist
    // only on rnaseq.genome (ribo-seq is single-end); secondary records in PE use
    // -f 320 (secondary + first-in-pair) so they are counted once per fragment.
    static final Map DEFAULT_COUNT_ARGS = [
        primary:   '-F 2304',          // not secondary, not supplementary
        secondary: '-f 256',
        unique:    '-q 255 -F 2304',   // STAR unique MAPQ, primary records only
    ].asImmutable()
    static final Map DEFAULT_PE_COUNT_ARGS = [
        paired_end_fragment:  '-f 64',
        paired_end_secondary: '-f 320',
    ].asImmutable()

    static Map route_count_args(Map params, String route) {
        def m = new LinkedHashMap(DEFAULT_COUNT_ARGS)
        if (route == 'rnaseq_genome') m.putAll(DEFAULT_PE_COUNT_ARGS)
        def user = route_params(params, route).samtools_count_arguments
        if (user instanceof Map) user.each { k, v -> m[k.toString()] = (v == null) ? '' : v.toString() }
        return m
    }

    // The MAPQ threshold implied by a counter's flags (BED-based counters need the
    // number itself). No -q ⇒ 0 (every record qualifies).
    static int count_mapq_threshold(String args) {
        def toks = shell_split(args ?: '')
        for (int i = 0; i < toks.size(); i++) {
            def t = toks[i]
            if (t in SAMTOOLS_MAPQ_OPTS && i + 1 < toks.size()) return toks[i + 1] as int
            def eq = t.indexOf('=')
            if (eq > 0 && t.substring(0, eq) in SAMTOOLS_MAPQ_OPTS) return t.substring(eq + 1) as int
        }
        return 0
    }

    // Shell block running the stats counters over `bam` concurrently (independent
    // reads of one file; wall time ≈ one pass). `bam` / `out_prefix` may contain
    // shell variables (e.g. "$l") — the text is inserted verbatim into the script.
    static String samtools_count_block(Map ca, boolean is_pe, String bam, String out_prefix, int cpus,
                                       boolean primary_secondary, boolean unique) {
        def frag = is_pe ? shell_quote_args((ca.paired_end_fragment ?: '').toString()) : ''
        def sec  = shell_quote_args(((is_pe ? ca.paired_end_secondary : null) ?: ca.secondary ?: '').toString())
        def prim = shell_quote_args((ca.primary ?: '').toString())
        def uniq = shell_quote_args((ca.unique ?: '').toString())
        def cmds = ['samtools view -c ' + frag + ' ' + bam + ' > ' + out_prefix + '.total.count']
        if (primary_secondary) {
            cmds << ('samtools view -c ' + prim + ' ' + frag + ' ' + bam + ' > ' + out_prefix + '.primary.count')
            cmds << ('samtools view -c ' + sec  + ' '        + bam + ' > ' + out_prefix + '.secondary.count')
        }
        if (unique) {
            cmds << ('samtools view -c ' + uniq + ' ' + frag + ' ' + bam + ' > ' + out_prefix + '.unique.count')
        }
        int t = Math.max(1, (cpus as int).intdiv(cmds.size()))
        return parallel_block(cmds.collect { it.replace('samtools view -c', 'samtools view -@ ' + t + ' -c') })
    }

    // Run shell commands concurrently and fail if any of them fails (a bare `wait`
    // always returns 0). The `${pids[@]+...}` form keeps `set -u` happy on bash 3.2.
    static String parallel_block(List<String> cmds) {
        def out = ['pids=()']
        cmds.each { c -> out << (c + ' & pids+=($!)') }
        out << 'for p in ${pids[@]+"${pids[@]}"}; do wait "$p"; done'
        return out.join('\n    ')
    }

    // ── validation ──────────────────────────────────────────────────────────
    static final List<String> SAMTOOLS_FILTER_OPTS_WITH_VALUE = [
        '-f', '--require-flags',
        '-F', '--excl-flags', '--exclude-flags',
        '-G',
        '-e', '--expr',
        '-d', '--tag',
        '-D', '--tag-file',
        '-L', '--target-file', '--region-file',
        '-r', '--read-group',
        '-R', '--read-group-file',
        '-N', '--qname-file',
        '-s', '--subsample', '--subsample-seed',
        '-m', '--min-qlen',
    ]
    static final List<String> SAMTOOLS_FILTER_FLAGS_NO_VALUE = ['-M', '--use-index', '-P', '--fetch-pairs']
    static final List<String> SAMTOOLS_MAPQ_OPTS = ['-q', '--min-MQ', '--min-mq']

    // Split a user argument string into tokens, honouring single and double quotes,
    // so a quoted filter expression (-e 'qlen>=25') survives as ONE token.
    static List<String> shell_split(String s) {
        def toks = []
        def cur = new StringBuilder()
        boolean started = false
        int i = 0
        int n = s.length()
        while (i < n) {
            String c = s.substring(i, i + 1)
            if (c == ' ' || c == '\t' || c == '\n' || c == '\r') {
                if (started) { toks << cur.toString(); cur = new StringBuilder(); started = false }
                i++
            }
            else if (c == "'" || c == '"') {
                int j = s.indexOf(c, i + 1)
                if (j < 0) throw new IllegalArgumentException("unterminated ${c} quote")
                cur.append(s.substring(i + 1, j)); started = true; i = j + 1
            }
            else { cur.append(c); started = true; i++ }
        }
        if (started) toks << cur.toString()
        return toks
    }

    // Re-quote every token for safe interpolation into a process script, so nothing
    // a user writes can reach the shell as syntax.
    static String shell_quote_args(String s) {
        if (s == null || s.trim().isEmpty()) return ''
        return shell_split(s).collect { "'" + it.replace("'", "'\\''") + "'" }.join(' ')
    }

    // Human-readable problems with one `samtools view` argument string; empty = valid.
    // `key` names the YAML key in messages. -q is rejected unless allow_mapq (the
    // quality filter owns it via mapping_quality_cutoff / unique_only; the `unique`
    // counter legitimately carries one).
    static List<String> validate_samtools_args(String args, String key, boolean allow_mapq = false) {
        def errors = []
        if (args == null) return errors
        if (args.toString() == 'true' || args.toString() == 'false') {
            return ["${key}: resolved to the boolean `${args}`. A bare `--${key} ''` on the " +
                    'command line is read by Nextflow as a flag, not an empty string — set the key ' +
                    'in your params YAML instead if you meant "no extra filtering".'].collect { it.toString() }
        }
        List<String> toks
        try { toks = shell_split(args.toString()) }
        catch (IllegalArgumentException e) { return ["${key}: ${e.message} in ${args.inspect()}".toString()] }
        boolean expect_value = false
        toks.each { t ->
            if (expect_value) { expect_value = false; return }
            if (!t.startsWith('-')) {
                errors << ("${key}: unexpected argument `${t}` — only options are allowed; RiboFlow " +
                           'supplies the input BAM and any region arguments itself.').toString()
                return
            }
            String name = t.contains('=') ? t.substring(0, t.indexOf('=')) : t
            if (SAMTOOLS_MAPQ_OPTS.contains(name)) {
                if (allow_mapq) { if (!t.contains('=')) expect_value = true }
                else errors << ("${key}: `${name}` is not allowed here — set `mapping_quality_cutoff` / " +
                                '`unique_only` instead, so the filter and the stats mode cannot disagree.').toString()
            }
            else if (SAMTOOLS_FILTER_OPTS_WITH_VALUE.contains(name)) {
                if (!t.contains('=')) expect_value = true
            }
            else if (!SAMTOOLS_FILTER_FLAGS_NO_VALUE.contains(name)) {
                errors << ("${key}: `${name}` is not an accepted read-filtering option. RiboFlow owns the " +
                           'output stream (-b/-h/-o/-@/-c/…); only record-selection options are allowed: ' +
                           (SAMTOOLS_FILTER_OPTS_WITH_VALUE + SAMTOOLS_FILTER_FLAGS_NO_VALUE).join(' ')).toString()
            }
        }
        if (expect_value) errors << "${key}: `${toks[-1]}` expects a value but none follows.".toString()
        return errors
    }

    static List<String> validate_count_args(Map params, String route) {
        def errors = []
        def user = route_params(params, route).samtools_count_arguments
        if (user == null) return errors
        def key = "${route}.samtools_count_arguments"
        def allowed = DEFAULT_COUNT_ARGS.keySet().toList() +
                      (route == 'rnaseq_genome' ? DEFAULT_PE_COUNT_ARGS.keySet().toList() : [])
        if (!(user instanceof Map)) {
            return ["${key}: must be a map with keys ${allowed.join(', ')}".toString()]
        }
        user.each { k, v ->
            if (!(k.toString() in allowed)) {
                def hint = k.toString().startsWith('paired_end')
                    ? ' (paired-end keys exist only under rnaseq.genome — ribo-seq is single-end)' : ''
                errors << "${key}.${k}: unknown key; allowed: ${allowed.join(', ')}${hint}".toString()
            } else {
                errors.addAll(validate_samtools_args(v?.toString(), "${key}.${k}".toString(), true))
            }
        }
        return errors
    }

    // All route-level problems, checked for every route regardless of which are
    // enabled, so a typo in a disabled block still surfaces immediately.
    static List<String> validate_routes(Map params) {
        def errors = legacy_filter_flag_errors(params)
        ROUTES.each { route ->
            errors.addAll(validate_samtools_args(route_filter_args(params, route),
                                                 "${route}.samtools_filter_arguments".toString()))
            errors.addAll(validate_count_args(params, route))
            def uo = route_params(params, route).unique_only
            if (uo != null) {
                try { as_bool(uo, true) }
                catch (IllegalArgumentException e) { errors << "${route}.unique_only: ${e.message}".toString() }
            }
            if (uo != null && !is_genome_route(route)) {
                errors << "${route}.unique_only: only the genome routes take this key (bowtie2 MAPQ is a confidence score; use mapping_quality_cutoff).".toString()
            }
        }
        return errors
    }

    // Non-fatal notes about the quality-filter configuration.
    static List<String> route_notes(Map params, List<String> active_routes) {
        def notes = []
        active_routes.findAll { is_genome_route(it) }.each { route ->
            def rp = route_params(params, route)
            if (rp.unique_only == null) {
                notes << ("${route}: `unique_only` is not set — inferred ${route_unique_only(params, route) ? 'unique-only' : 'multi-mapper'} " +
                          "mode from mapping_quality_cutoff=${route_mapq_cutoff(params, route)}. Set `${route}.unique_only: true|false` explicitly.").toString()
            }
            else if (route_unique_only(params, route) && rp.mapping_quality_cutoff != null &&
                     route_mapq_cutoff(params, route) != STAR_UNIQUE_MAPQ) {
                notes << ("${route}: unique_only is true, so `samtools view -q ${STAR_UNIQUE_MAPQ}` is applied and " +
                          "mapping_quality_cutoff=${route_mapq_cutoff(params, route)} is ignored.").toString()
            }
        }
        return notes
    }

    // The removed integer keys. Fail loudly with the translation rather than letting a
    // stale YAML silently fall through to the default mask.
    static List<String> legacy_filter_flag_errors(Map params) {
        def errors = []
        def check = { Map holder, String old_key, String old_path, String new_path ->
            if (holder != null && holder.containsKey(old_key)) {
                errors << ("`${old_path}` was removed. Use `${new_path}: \"-F ${holder[old_key]}\"` " +
                           'instead — it takes any samtools view read-filtering arguments, not just a ' +
                           'single -F mask.').toString()
            }
        }
        def rnaseq = (params.rnaseq ?: [:]) as Map
        check((params.genome ?: [:]) as Map, 'ribo_filter_flags',
              'genome.ribo_filter_flags', 'genome.samtools_filter_arguments')
        check((rnaseq.genome ?: [:]) as Map, 'filter_flags',
              'rnaseq.genome.filter_flags', 'rnaseq.genome.samtools_filter_arguments')
        check(rnaseq, 'filter_flags',
              'rnaseq.filter_flags', 'rnaseq.genome.samtools_filter_arguments')
        return errors
    }

    // ── per-sample → per-lane fan-out ───────────────────────────────────────
    // Processes that split one merged sample file into per-lane files emit them as
    // a list; these map them back onto the lane metas by the `<id>.<lane>.` prefix.
    static List as_list(Object x) {
        return x == null ? [] : (x instanceof List ? x : [x])
    }

    static Object lane_file(Object files, Map meta) {
        def pfx = "${meta.id}.${meta.lane}.".toString()
        def hit = as_list(files).find { f -> f.toString().tokenize('/')[-1].startsWith(pfx) }
        if (hit == null) throw new IllegalStateException("no per-lane file for ${pfx}* among ${as_list(files)*.toString()}")
        return hit
    }

    static List lane_pairs(Object files, List metas) {
        return metas.collect { m -> [m, lane_file(files, m)] }
    }

    // Several per-sample file lists → one [meta, f1, f2, ...] tuple per lane.
    static List lane_tuples(List file_lists, List metas) {
        return metas.collect { m -> [m] + file_lists.collect { fl -> lane_file(fl, m) } }
    }

    static List lane_ids(List metas) {
        return metas.collect { m -> "${m.id}.${m.lane}".toString() }
    }
}
