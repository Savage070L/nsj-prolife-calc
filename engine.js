(function () {
    var _z = "eyJtIjpbMCwwLDAsMCwwLDAsMCwwLDAsMCwwLDAsMCwwLDAsMCwxLjM5LDEuNzUsMS44NCwxLjg0LDEuODUsMS44NywxLjg4LDEuODgsMS44OSwxLjkxLDEuOTIsMS45MywxLjkzLDEuOTQsMS45NiwxLjk3LDEuOTcsMS45OCwyLDIuMTEsMi4yNywyLjQ0LDIuNjYsMi44OSwzLjE0LDMuMzksMy42NiwzLjk3LDQuMyw0LjY5LDUuMTQsNS43Miw2LjM0LDcuMDMsNy43NCw4LjQ2LDkuMjEsMTAsMTAuODgsMTEuODgsMTMuMDMsMTQuMzUsMTUuODEsMTcuNDEsMTkuMDcsMjAuODMsMjIuNzIsMjQuOCwyNy4xNCwyOS44LDMyLjgyLDM2LjIyLDQwLjA1LDQ0LjMzLDQ5LjExLDU0LjI3LDU5Ljk3LDY2LjI3LDczLjIzLDgwLjkyLDg5LjQyLDk4LjgyLDEwOS4yLDEyMC42NywxMzMuMzUsMTQ3LjM2LDE2Mi44NCwxNzkuOTUsMTk4Ljg2LDIxOS43NiwyNDIuODUsMjY4LjM3LDI5Ni41NywzMjcuNzMsMzYyLjE3LDQwMC4yMiw0NDIuMjcsNDg4Ljc0LDU0MC4wOSw1OTYuODQsNjU5LjU1LDcyOC44NSw4MDUuNDMsODkwLjA2LDk4My41OF0sImYiOlswLDAsMCwwLDAsMCwwLDAsMCwwLDAsMCwwLDAsMCwwLDEuMDYsMS4wOCwwLjk4LDEuMDEsMS4wMiwxLjA0LDEuMDYsMS4wOCwxLjEsMS4xMSwxLjEzLDEuMTUsMS4xNywxLjIsMS4yNCwxLjI4LDEuMzMsMS40MSwxLjQ4LDEuNTUsMS42NiwxLjc1LDEuODUsMS45NiwyLjA2LDIuMTksMi4zNCwyLjUxLDIuNywyLjkxLDMuMTQsMy40LDMuNjgsNCw0LjMzLDQuNjksNS4wNiw1LjQ0LDUuODYsNi4yOSw2LjgsNy40Niw4LjIxLDkuMDQsOS45MywxMC44OSwxMS45LDEzLjA0LDE0LjMyLDE1Ljc5LDE3LjUsMTkuNDUsMjEuNzYsMjQuNDYsMjcuNjYsMzAuOTQsMzQuNjEsMzguNzIsNDMuMzIsNDguNDYsNTQuMjEsNjAuNjQsNjcuODQsNzUuODksODQuOSw5NC45OCwxMDYuMjUsMTE4Ljg2LDEzMi45NywxNDguNzUsMTY2LjQsMTg2LjE1LDIwOC4yNCwyMzIuOTUsMjYwLjYsMjkxLjUzLDMyNi4xMywzNjQuODMsNDA4LjEzLDQ1Ni41Nyw1MTAuNzYsNTcxLjM4LDYzOS4xOSw3MTUuMDUsNzk5LjkxXSwiY20iOlswLDAsMCwwLDAsMCwwLDAsMCwwLDAsMCwwLDAsMCwwLDAuMzIsMC4zMiwwLjMyLDAuMzIsMC4zMiwwLjM1LDAuNDMsMC40MywwLjQ1LDAuNDUsMC40NSwwLjQyLDAuNDUsMC40NSwwLjUzLDAuNjYsMC43NSwwLjkxLDEuMTIsMS4yOSwxLjU1LDEuODYsMi4yOCwyLjc3LDMuMzIsMy43NSw0LjM1LDQuNTgsNC44Myw1LjI3LDUuNzEsNi4xNCw2LjksNy43Miw4LjQyLDkuMzEsMTAuMzEsMTEuNTcsMTIuNTksMTMuNywxNS4wNCwxNi41OCwxNy42MSwxOC41NiwxOS4zOSwyMy4xNywyNy4wOSwyOC44MSwzMC41NywzMi40MSwzNS45NywzOS45Miw0NC4zLDQ5LjE2LDU0LjU2LDYwLjU1LDY3LjIsNzQuNTgsODIuNzcsOTEuODYsMTAxLjk0LDExMy4xMywxMjUuNTUsMTM5LjMzLDE1NC42MywxNzEuNjEsMTkwLjQ1LDIxMS4zNiwyMzQuNTYsMjYwLjMxLDI4OC44OSwzMjAuNiwzNTUuOCwzOTQuODYsNDM4LjIxLDQ4Ni4zMiw1MzkuNzEsNTk4Ljk2LDY2NC43MSw3MzcuNjgsODE4LjY2LDkwOC41MywxMDA4LjI3LDExMTguOTYsMTI0MS44XSwiY2YiOlswLDAsMCwwLDAsMCwwLDAsMCwwLDAsMCwwLDAsMCwwLDAuNDYsMC40NiwwLjQ2LDAuNDYsMC40NiwwLjQ4LDAuNjIsMC42OCwwLjc0LDAuODMsMC45MywwLjk5LDEuMTMsMS4yMiwxLjMsMS40MSwxLjYyLDEuNzMsMS45MSwyLjEsMi4zMywyLjU1LDIuNzcsMi45NywzLjIyLDMuNSwzLjc0LDMuODYsNC4xLDQuNDIsNC43Nyw1LjEyLDUuNTcsNS45OSw2LjM1LDYuNiw3LjE5LDcuNjksNy45Nyw4LjMyLDguODEsOS4yNyw5LjY0LDEwLjM5LDEwLjk5LDEyLjgzLDE0Ljc4LDE1Ljc5LDE2LjYsMTcuNSwxOS4yMywyMS4xMywyMy4yMSwyNS41LDI4LjAyLDMwLjc4LDMzLjgyLDM3LjE2LDQwLjgzLDQ0Ljg2LDQ5LjI5LDU0LjE1LDU5LjQ5LDY1LjM2LDcxLjgxLDc4Ljg5LDg2LjY3LDk1LjIyLDEwNC42MSwxMTQuOTMsMTI2LjI3LDEzOC43MywxNTIuNDIsMTY3LjQ2LDE4My45OCwyMDIuMTMsMjIyLjA3LDI0My45OCwyNjguMDUsMjk0LjQ5LDMyMy41NCwzNTUuNDYsMzkwLjUzLDQyOS4wNiw0NzEuMzldfQ==";
    var _q0 = JSON.parse(atob(_z));
    var _h0 = 0.07, _h1 = 0.12, _w0 = 100, _w2 = 0.01;
    var _fa = { annual: 1.0, semiannual: 0.51, quarterly: 0.2575, monthly: 0.0875 };
    var _ex = { "3": { K: 0.11923, L: 0.02, M: 0.021, N: 0.02, O: 0.0 }, "4": { K: 0.1523, L: 0.02, M: 0.025, N: 0.02, O: 0.0 }, "5": { K: 0.18538, L: 0.04, M: 0.027, N: 0.025, O: 0.0 }, "6": { K: 0.2946, L: 0.1434, M: 0.029, N: 0.03, O: 0.003 }, "7": { K: 0.3387, L: 0.1623, M: 0.030, N: 0.03, O: 0.003 }, "8": { K: 0.3828, L: 0.1812, M: 0.03, N: 0.03, O: 0.003 }, "9": { K: 0.4269, L: 0.2001, M: 0.03, N: 0.03, O: 0.003 }, "10": { K: 0.501, L: 0.249, M: 0.03, N: 0.04, O: 0.003 }, "11": { K: 0.5451, L: 0.2679, M: 0.03, N: 0.04, O: 0.003 }, "12": { K: 0.5892, L: 0.2868, M: 0.03, N: 0.04, O: 0.003 }, "13": { K: 0.6333, L: 0.3057, M: 0.03, N: 0.04, O: 0.003 }, "14": { K: 0.6774, L: 0.3246, M: 0.03, N: 0.04, O: 0.003 }, "15": { K: 0.7215, L: 0.3435, M: 0.025, N: 0.04, O: 0.003 } };
    var _rd = {
        accidental_death: { t: 0.001008, e: 0.25, q: 0.4 },
        traffic_death: { t: 0.000336, e: 0.25, q: 0.4 },
        disability_accident_lumpsum: { t: 0.00048, e: 0.25, q: 0.4 },
        disability_any_lumpsum: { t: 0.001392, e: 0.25, q: 0.4 },
        trauma: { t: 0.0045, e: 0.05, q: 0.25 },
        trauma_extra: { t: 0.005643, e: 0.05, q: 0.25 },
        temporary_disability: { t: 0.002, e: 0.05, q: 0.25 },
        hospitalization: { t: 0.000786, e: 0.05, q: 0.25 }
    };
    /* Build mortality commutation table */
    function _b1(r, i) {
        var v = 1 / (1 + i), u = [], l = [], D = [], N = [], C = [], M = [], x;
        for (x = 0; x < 101; x++)u[x] = r[x] / 1000;
        u[100] = 1.0;
        l[0] = 1000000;
        for (x = 0; x < 100; x++)l[x + 1] = l[x] * (1 - u[x]);
        for (x = 0; x < 101; x++)D[x] = l[x] * Math.pow(v, x);
        N[100] = D[100];
        for (x = 99; x >= 0; x--)N[x] = D[x] + N[x + 1];
        for (x = 0; x < 100; x++)C[x] = (l[x] - l[x + 1]) * Math.pow(v, x + 1);
        C[100] = (l[100] >= 1) ? 1 : 0;
        M[100] = C[100];
        for (x = 99; x >= 0; x--)M[x] = C[x] + M[x + 1];
        return { D: D, N: N, C: C, M: M, l: l, q: u };
    }

    /* Build CI commutation table (double decrement: death + CI) */
    function _b2(rd, rc, i) {
        var v = 1 / (1 + i), qd = [], qc = [], x;
        for (x = 0; x < 101; x++) { qd[x] = rd[x] / 1000; qc[x] = rc[x] / 1000; }
        qd[100] = 1.0;
        for (x = 0; x < 101; x++) {
            var s = qd[x] + qc[x];
            if (s > 1) { qd[x] = qd[x] / s; qc[x] = qc[x] / s; }
        }
        var l = [], D = [], N = [], AG = [];
        l[0] = 1000000;
        for (x = 0; x < 100; x++)l[x + 1] = Math.max(l[x] * (1 - qd[x] - qc[x]), 0);
        for (x = 0; x < 101; x++)D[x] = l[x] * Math.pow(v, x);
        N[100] = D[100];
        for (x = 99; x >= 0; x--)N[x] = D[x] + N[x + 1];
        for (x = 0; x < 101; x++)AG[x] = D[x] * qc[x];
        return { D: D, N: N, l: l, qc: qc };
    }

    /* Clamp expense key */
    function _ek(n) { return String(Math.min(Math.max(n, 3), 15)); }

    /* Get expense params G2,G3,G6,G7 */
    function _gx(n, t) {
        var g6 = _ex[_ek(n)].N, g7 = _ex[_ek(n)].O;
        var g10 = Math.ceil(_ex[_ek(n)].M * 100) / 100;
        var g2, g3;
        if (t === 1) { g2 = g10; g3 = g10; }
        else { var kt = _ek(t); g2 = Math.ceil(_ex[kt].K * 100) / 100; g3 = Math.ceil(_ex[kt].L * 100) / 100; }
        return { g2: g2, g3: g3, g6: g6, g7: g7 };
    }

    /* Main premium calculation: actuarial values + BP */
    function _p1(tb, x, n, t) {
        var Dx = tb.D[x], Dxn = tb.D[x + n], Dx1 = tb.D[x + 1];
        var Mx = tb.M[x], Mxn = tb.M[x + n];
        var Nx = tb.N[x], Nxn = tb.N[x + n], Nxt = tb.N[x + t];
        var Ax = (Mx - Mxn + Dxn) / Dx;
        var an = (Nx - Nxn) / Dx;
        var at = (Nx - Nxt) / Dx;
        var NP = (at > 0) ? Ax / at : 0;
        var ge = _gx(n, t);
        var num = Ax + ge.g7 * an;
        var den = at - ge.g6 * at - (ge.g2 + ge.g3 * Dx1 / Dx);
        var BP = (den > 0) ? num / den : 0;
        return { Ax: Ax, an: an, at: at, BP: BP, NP: NP, ge: ge, Dx: Dx, Dx1: Dx1, Dxn: Dxn, Mx: Mx, Mxn: Mxn, Nx: Nx, Nxn: Nxn, Nxt: Nxt };
    }

    /* Reserves table */
    function _r1(tb, x, n, t, BP, ge, SA) {
        var rv = [], k;
        for (k = 1; k <= n; k++) {
            var Dxk = tb.D[x + k], Dxn = tb.D[x + n];
            var Mxk = tb.M[x + k], Mxn = tb.M[x + n];
            var Nxk = tb.N[x + k], Nxn = tb.N[x + n], Nxt = tb.N[x + t];
            if (Dxk === 0) { rv.push({ year: k, reserve: 0, surrender: 0, reduced_sa: 0 }); continue; }
            var Axk = (Mxk - Mxn + Dxn) / Dxk;
            var ank = (Nxk - Nxn) / Dxk;
            var atk = (k < t) ? (Nxk - Nxt) / Dxk : 0;
            var al = (k === 1) ? ge.g3 : 0;
            var rr = Axk + ge.g7 * ank - BP * (atk - ge.g6 * atk - al);
            var sr = rr - (1 - rr) * _w2;
            var res = rr * SA;
            var sur = Math.max(sr * SA, 0);
            if (k === n) sur = SA;
            var reduced_sa = (Axk > 0) ? Math.round(sur) / Axk : 0;
            rv.push({ year: k, reserve: Math.round(res), surrender: Math.round(sur), reduced_sa: Math.round(reduced_sa) });
        }
        return rv;
    }

    /* Simple rider premium */
    function _p2(rk, rs, n, ff, sg) {
        var rc = _rd[rk];
        if (!rc) return null;
        var bt = rc.t;
        var gt = Math.round((bt * (1 + rc.e)) / (1 - rc.q) * 10000) / 10000;
        if (sg) return { gross_tariff: gt, rider_sum: rs, rider_premium: Math.round(gt * rs * n) };
        return { gross_tariff: gt, rider_sum: rs, rider_premium: Math.round(gt * rs * ff) };
    }

    /* CI rider premium (actuarial, double decrement) */
    function _p3(qx, qci, i, x, n, t, ge, ff, cs, sg) {
        var ctb = _b2(qx, qci, i);
        var Dx0 = ctb.D[x];
        if (Dx0 === 0) return { BP_ci: 0, ci_sum: cs, premium: 0 };
        /* A_ci = sum D(x+k)/D(x) * q_ci(x+k) for k=0..n-1 */
        var A_ci = 0, k;
        for (k = 0; k < n; k++) {
            var age_k = x + k;
            if (age_k > _w0) break;
            A_ci += ctb.D[age_k] / Dx0 * ctb.qc[age_k];
        }
        var Nxci = ctb.N[x], Nxnci = ctb.N[x + n], Nxtci = ctb.N[x + t];
        var Dx1ci = ctb.D[x + 1];
        var ax_n_ci = (Nxci - Nxnci) / Dx0;
        var ax_t_ci = (Nxci - Nxtci) / Dx0;
        var num = A_ci + ge.g7 * ax_n_ci;
        var den = ax_t_ci - ge.g6 * ax_t_ci - (ge.g2 + ge.g3 * Dx1ci / Dx0);
        var BP_ci = (den > 0) ? num / den : 0;
        BP_ci = Math.round(BP_ci * 10000) / 10000;
        var pm;
        if (sg) pm = Math.round(BP_ci * cs);
        else pm = Math.round(BP_ci * cs * ff);
        return { BP_ci: BP_ci, A_ci: A_ci, ax_n_ci: ax_n_ci, ax_t_ci: ax_t_ci, ci_sum: cs, premium: pm };
    }

    /* Premium waiver: PV of future premiums × disability_accident tariff */
    function _p4w(annual_premium_main, ff, freq, sg, n) {
        if (sg || freq === 'single') return 0;
        var ap = Math.round(annual_premium_main);
        var rc = _rd.disability_accident_lumpsum;
        var j6 = Math.round((rc.t * (1 + rc.e) / (1 - rc.q)) * 10000) / 10000;
        var r_sum = 0;
        for (var k = 0; k < n - 1; k++) {
            var q_k = (n - 1 - k) * ap;
            r_sum += Math.round(q_k * j6);
        }
        if (n <= 1) return 0;
        return Math.round(r_sum / (n - 1) * ff);
    }

    /* Helper: calculate one rider premium for iteration */
    function _rp1(rk, rs, x, n, t, gn, freq, ff, sg, annual_pm, qx, qci, ge) {
        if (rk === 'premium_waiver') {
            return _p4w(annual_pm, ff, freq, sg, n);
        } else if (rk === 'critical_illness') {
            var ci = _p3(qx, qci, _h0, x, n, t, ge, ff, rs, sg);
            return ci.premium;
        } else {
            var sr = _p2(rk, rs, n, ff, sg);
            return sr ? sr.rider_premium : 0;
        }
    }

    /* Main calculate function */
    function _ca(p) {
        var db = new Date(p.dob), td = new Date();
        var ag = td.getFullYear() - db.getFullYear();
        var m = td.getMonth() - db.getMonth();
        if (m < 0 || (m === 0 && td.getDate() < db.getDate())) ag--;
        var x = ag, n = parseInt(p.term), fr = p.frequency || 'annual', md = p.mode || 'sa_to_premium';
        var sg = (fr === 'single'), i = (sg && n === 3) ? _h1 : _h0;
        var gn = (p.gender === 'male') ? 'm' : 'f';
        var qx = _q0[gn], ff = sg ? 1 : _fa[fr];
        var exit_age = x + n;

        /* Age validation */
        if (exit_age > 70) {
            return {
                success: false,
                error: '\u0412\u043e\u0437\u0440\u0430\u0441\u0442 \u0437\u0430\u0441\u0442\u0440\u0430\u0445\u043e\u0432\u0430\u043d\u043d\u043e\u0433\u043e \u043d\u0430 \u043c\u043e\u043c\u0435\u043d\u0442 \u043e\u043a\u043e\u043d\u0447\u0430\u043d\u0438\u044f \u0434\u043e\u0433\u043e\u0432\u043e\u0440\u0430 \u0441\u043e\u0441\u0442\u0430\u0432\u0438\u0442 ' + exit_age + ' \u043b\u0435\u0442. \u041c\u0430\u043a\u0441\u0438\u043c\u0430\u043b\u044c\u043d\u043e \u0434\u043e\u043f\u0443\u0441\u0442\u0438\u043c\u044b\u0439 \u0432\u043e\u0437\u0440\u0430\u0441\u0442 \u043d\u0430 \u0432\u044b\u0445\u043e\u0434\u0435 \u2014 70 \u043b\u0435\u0442. \u0423\u043c\u0435\u043d\u044c\u0448\u0438\u0442\u0435 \u0441\u0440\u043e\u043a \u0441\u0442\u0440\u0430\u0445\u043e\u0432\u0430\u043d\u0438\u044f \u0438\u043b\u0438 \u043f\u0440\u043e\u0432\u0435\u0440\u044c\u0442\u0435 \u0434\u0430\u0442\u0443 \u0440\u043e\u0436\u0434\u0435\u043d\u0438\u044f.',
                validation: 'exit_age_exceeded',
                age: x,
                exit_age: exit_age
            };
        }

        var tb = _b1(qx, i);
        var t = sg ? 1 : n;
        var r = _p1(tb, x, n, t);
        var BP = r.BP, ge = r.ge, NP = r.NP;

        /* Parse rider groups */
        var ri_in = p.riders || {};
        var g1k = ri_in.group1 || null;
        var g2k = ri_in.group2 || null;
        var g3 = ri_in.group3 || {};
        var sa_linked = [];
        if (g1k) sa_linked.push(g1k);
        if (g2k) sa_linked.push(g2k);
        var fixed_r = [];
        var _g3keys = ['trauma', 'trauma_extra', 'temporary_disability', 'hospitalization', 'critical_illness'];
        for (var gi = 0; gi < _g3keys.length; gi++) {
            var _gk = _g3keys[gi];
            var _gs = g3[_gk] || {};
            if (_gs.enabled) fixed_r.push([_gk, parseFloat(_gs.sum || 0)]);
        }

        /* CI data for rider calcs */
        var gci = (p.gender === 'male') ? 'cm' : 'cf';
        var qci_raw = _q0[gci];

        var SA, ap, gp;

        if (md === 'premium_to_sa') {
            /* === Premium -> SA (iterative) === */
            var total_pm_input = parseFloat(p.premium);

            /* Step 1: fixed rider premiums (don't depend on SS) */
            var fixed_total = 0;
            for (var fi = 0; fi < fixed_r.length; fi++) {
                var fk = fixed_r[fi][0], fs = fixed_r[fi][1];
                var fp = _rp1(fk, fs, x, n, t, gn, fr, ff, sg, 0, qx, qci_raw, ge);
                fixed_total += fp;
            }

            /* Step 2: initial SS guess */
            var remaining = total_pm_input - fixed_total;
            if (sg) SA = (BP > 0) ? remaining / BP : 0;
            else SA = (BP > 0) ? remaining / (BP * ff) : 0;

            /* Step 3: iterate */
            for (var it = 0; it < 100; it++) {
                var apm_iter = Math.round(BP * SA);
                var sa_total = 0;
                for (var si = 0; si < sa_linked.length; si++) {
                    var sk = sa_linked[si];
                    var rp_iter = _rp1(sk, SA, x, n, t, gn, fr, ff, sg, apm_iter, qx, qci_raw, ge);
                    sa_total += rp_iter;
                }
                remaining = total_pm_input - fixed_total - sa_total;
                if (remaining <= 0) remaining = 0;
                var new_sa;
                if (sg) new_sa = (BP > 0) ? remaining / BP : 0;
                else new_sa = (BP > 0) ? remaining / (BP * ff) : 0;
                if (Math.abs(new_sa - SA) < 1) { SA = new_sa; break; }
                SA = new_sa;
            }
            SA = Math.round(SA);
            ap = Math.round(BP * SA);
            gp = sg ? ap : Math.round(BP * SA * ff);
        } else {
            /* === SA -> Premium === */
            SA = parseFloat(p.sum_assured);
            ap = Math.round(BP * SA);
            gp = sg ? ap : Math.round(BP * SA * ff);
        }

        /* Reserves */
        var rv = _r1(tb, x, n, t, BP, ge, SA);

        /* Calculate riders */
        var ri = {}, rt = 0;

        /* SA-linked riders (group1 + group2) */
        for (var ai = 0; ai < sa_linked.length; ai++) {
            var ak = sa_linked[ai];
            if (ak === 'premium_waiver') {
                var wp = _p4w(ap, ff, fr, sg, n);
                ri.premium_waiver = { waiver_premium: wp, annual_payment: Math.round(ap) };
                rt += wp;
            } else {
                var ar = _p2(ak, SA, n, ff, sg);
                if (ar) { ri[ak] = ar; rt += ar.rider_premium; }
            }
        }

        /* Fixed riders (group3) */
        for (var bi = 0; bi < fixed_r.length; bi++) {
            var bk = fixed_r[bi][0], bs = fixed_r[bi][1];
            if (bk === 'critical_illness') {
                var ci = _p3(qx, qci_raw, i, x, n, t, ge, ff, bs, sg);
                ri.critical_illness = ci;
                rt += ci.premium;
            } else {
                var br = _p2(bk, bs, n, ff, sg);
                if (br) { ri[bk] = br; rt += br.rider_premium; }
            }
        }

        var net_premium = NP * SA;

        return {
            success: true,
            age: x,
            exit_age: exit_age,
            term: n,
            payment_term: t,
            frequency: fr,
            mode: md,
            sum_assured: SA,
            BP: Math.round(BP * 1000000) / 1000000,
            NP_rate: NP,
            BP_rate: BP,
            annual_premium: ap,
            gross_premium: gp,
            net_premium: net_premium,
            interest_rate: i,
            freq_factor: ff,
            reserves: rv,
            riders: ri,
            riders_total: rt,
            total_premium: gp + rt,
            sa_linked_riders: sa_linked,
            fixed_riders: fixed_r.map(function (f) { return f[0]; }),
            Ax_n: r.Ax,
            ax_n: r.an,
            ax_t: r.at,
            Dx: r.Dx,
            Dx1: r.Dx1,
            Dxn: r.Dxn,
            Mx: r.Mx,
            Mxn: r.Mxn,
            Nx: r.Nx,
            Nxn: r.Nxn,
            Nxt: r.Nxt,
            G2: ge.g2,
            G3: ge.g3,
            G6: ge.g6,
            G7: ge.g7
        };
    }

    window.NSJEngine = { calculate: _ca };
})();
