use crate::models::short_weierstrass_jacobian::GroupAffine;
use ark_ff::{fields::BitIteratorBE, Field, PrimeField, SquareRootField, Zero};

pub mod bls12;
pub mod bn;
pub mod bw6;
pub mod mnt4;
pub mod mnt6;
pub mod short_weierstrass_jacobian;
pub mod twisted_edwards_extended;

/// Model composed by a base and scalar field definitions
pub trait ModelParameters: Send + Sync + 'static {
    /// Base field of the model with efficient square-root implementation.
    type BaseField: Field + SquareRootField;
    /// Finite scalar field of the model
    type ScalarField: PrimeField + SquareRootField + Into<<Self::ScalarField as PrimeField>::BigInt>;
}

/// Model defined as short-weierstrass form
///
/// # Definition
///
/// `y^2 = x^3 + a · x + b`
pub trait SWModelParameters: ModelParameters {
    /// Coefficient `a` of the weierstrass equation
    const COEFF_A: Self::BaseField;
    /// Coefficient `b` of the weierstrass equation
    const COEFF_B: Self::BaseField;
    /// Cofactor that defines the subgroup of the curve
    const COFACTOR: &'static [u64];
    /// Inverse of the cofactor
    const COFACTOR_INV: Self::ScalarField;
    /// Coefficients of the base point of the curve
    const AFFINE_GENERATOR_COEFFS: (Self::BaseField, Self::BaseField);

    /// Multiply the provided field element by the weierstrass coefficient `a`
    #[inline(always)]
    fn mul_by_a(elem: &Self::BaseField) -> Self::BaseField {
        let mut copy = *elem;
        copy *= &Self::COEFF_A;
        copy
    }

    /// Sum the provided field element with the weierstrass coefficient `b`
    #[inline(always)]
    fn add_b(elem: &Self::BaseField) -> Self::BaseField {
        if !Self::COEFF_B.is_zero() {
            let mut copy = *elem;
            copy += &Self::COEFF_B;
            return copy;
        }
        *elem
    }

    /// Evaluate if the provided affine point is in the correct subgroup
    /// relative to the cofactor of the model.
    fn is_in_correct_subgroup_assuming_on_curve(item: &GroupAffine<Self>) -> bool
    where
        Self: Sized,
    {
        item.mul_bits(BitIteratorBE::new(Self::ScalarField::characteristic()))
            .is_zero()
    }
}

/// Model defined as twisted-edwards form
///
/// # Definition
///
/// `a · x^2 + y^2 = 1 + d · x^2 · y^2`
pub trait TEModelParameters: ModelParameters {
    /// Coefficient `a` of the twisted-edwards equation
    const COEFF_A: Self::BaseField;
    /// Coefficient `d` of the twisted-edwards equation
    const COEFF_D: Self::BaseField;
    /// Cofactor that defines the subgroup of the curve
    const COFACTOR: &'static [u64];
    /// Inverse of the cofactor
    const COFACTOR_INV: Self::ScalarField;
    /// Coefficients of the base point of the curve
    const AFFINE_GENERATOR_COEFFS: (Self::BaseField, Self::BaseField);

    /// Montgomery model with birational equivalence to this model
    type MontgomeryModelParameters: MontgomeryModelParameters<BaseField = Self::BaseField>;

    /// Multiply the provided field element by the twisted-edwards coefficient
    /// `a`
    #[inline(always)]
    fn mul_by_a(elem: &Self::BaseField) -> Self::BaseField {
        let mut copy = *elem;
        copy *= &Self::COEFF_A;
        copy
    }
}

/// Model defined as montgomery form
///
/// # Definition
///
/// `b · y^2 = x^3 + a · x^2 + x`
pub trait MontgomeryModelParameters: ModelParameters {
    /// Coefficient `a` of the montgomery equation
    const COEFF_A: Self::BaseField;
    /// Coefficient `b` of the montgomery equation
    const COEFF_B: Self::BaseField;

    /// Twisted-Edwards model with birational equivalence to this model
    type TEModelParameters: TEModelParameters<BaseField = Self::BaseField>;
}
