use libm;

// largest argument for exp, which means that exp(x) = NaN
pub(crate) const MAXLOG: f64 = 709.782712893383973096206318587;
// the machine roundoff error
pub(crate) const MACHEP: f64 = 1.11022302462515654042363166809e-16;
pub(crate) const BIG: f64 = 4.503599627370496e15;
pub(crate) const BIGINV: f64 = 2.22044604925031308085e-16;
pub(crate) const SQRT2: f64 = 1.41421356237309504880;

// The mathmatical functions.
// Reference:
// [nr] Numerical recipes in C (3nd ed) the art of scientific computing
// [AS] Milton Abramowitz, Irene A. Stegun Handbook of Mathematical Functions with Formulas, Graphs, and Mathematical Tables

#[inline(always)]
pub(crate) fn ln(x: f64) -> f64 {
    x.ln()
}

#[inline(always)]
pub(crate) fn exp(x: f64) -> f64 {
    x.exp()
}

#[inline(always)]
/// return x^y
pub(crate) fn powf(x: f64, y: f64) -> f64 {
    x.powf(y)
}

#[inline(always)]
/// return x^y
pub(crate) fn powi(x: f64, y: i32) -> f64 {
    x.powi(y)
}

#[inline(always)]
pub(crate) fn sqrt(x: f64) -> f64 {
    x.sqrt()
}

#[inline(always)]
pub(crate) fn erfc(x: f64) -> f64 {
    libm::erfc(x)
}

#[inline(always)]
pub(crate) fn erf(x: f64) -> f64 {
    libm::erf(x)
}

#[inline(always)]
pub(crate) fn lgam(x: f64) -> f64 {
    x.ln_gamma().0
}

#[inline(always)]
pub(crate) fn abs(x: f64) -> f64 {
    x.abs()
}

// incomplete complementary gamma function
// https://root.cern.ch/doc/v614/SpecFuncCephes_8cxx_source.html
pub(crate) fn igamc(a: f64, x: f64) -> f64 {
    if a <= 0.0 {
        return 0.0;
    }
    if x <= 0.0 {
        return 1.0;
    }

    if (x < 1.0) || (x < a) {
        return 1.0 - igam(a, x);
    }

    let mut ax = a * ln(x) - x - lgam(a);

    if ax < -MAXLOG {
        return 0.0;
    }
    ax = exp(ax);

    /* continued fraction */
    let mut y = 1.0 - a;
    let mut z = x + y + 1.0;
    let mut c = 0.0;
    let mut pkm2 = 1.0;
    let mut qkm2 = x;
    let mut pkm1 = x + 1.0;
    let mut qkm1 = z * x;
    let mut ans = pkm1 / qkm1;
    let mut t;
    let mut r;
    loop {
        c += 1.0;
        y += 1.0;
        z += 2.0;
        let yc = y * c;
        let pk = pkm1 * z - pkm2 * yc;
        let qk = qkm1 * z - qkm2 * yc;
        if qk != 0.0 {
            r = pk / qk;
            t = abs((ans - r) / r);
            ans = r;
        } else {
            t = 1.0;
        }
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;
        if abs(pk) > BIG {
            pkm2 *= BIGINV;
            pkm1 *= BIGINV;
            qkm2 *= BIGINV;
            qkm1 *= BIGINV;
        }
        if t <= MACHEP {
            break;
        }
    }

    ans * ax
}

// incomplete gamma function
pub(crate) fn igam(a: f64, x: f64) -> f64 {
    // LM: for negative values returns 1.0 instead of zero
    // This is correct if a is a negative integer since Gamma(-n) = +/- inf
    if a <= 0.0 {
        return 1.0;
    }

    if x <= 0.0 {
        return 0.0;
    }

    if (x > 1.0) && (x > a) {
        return 1.0 - igamc(a, x);
    }

    /* Compute  x**a * exp(-x) / gamma(a)  */
    let mut ax = a * ln(x) - x - lgam(a);
    if ax < -MAXLOG {
        return 0.0;
    }

    ax = exp(ax);

    /* power series */
    let mut r = a;
    let mut c = 1.0;
    let mut ans = 1.0;

    loop {
        r += 1.0;
        c *= x / r;
        ans += c;
        if c / ans <= MACHEP {
            break;
        }
    }

    ans * ax / a
}

pub(crate) fn normal(x: f64) -> f64 {
    if x.is_sign_positive(){
        0.5 * (1.0 + erf(x / SQRT2))
    }else{
        0.5 * (1.0 - erf(-x / SQRT2))
    }
}
