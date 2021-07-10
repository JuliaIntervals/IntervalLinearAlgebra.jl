function solve(A, b,
               solver::AbstractIterativeSolver,
               precondition::AbstractPrecondition=_default_precondition(A, solver))

    A, b = precondition(A, b)
    x = enclose(A, b)

    solver(A, b, x)

    return x
end

function solve(A, b,
               solver::AbstractDirectSolver,
               precondition::AbstractPrecondition=_default_precondition(A, solver))

    A, b = precondition(A, b)
    return solver(A, b)
end


# fallback
function solve(A, b,
               solver::AbstractLinearSolver=_default_solver(),
               precondition::AbstractPrecondition=_default_precondition(A, solver))

    A, b = precondition(A, b)
    return solver(A, b)
end

## Default settings
_default_solver() = GaussianElimination()

function _default_precondition(A, ::GaussianElimination)
    if is_strictly_diagonally_dominant(A) || is_M_matrix(A)
        return NoPrecondition()
    else
        return InverseMidpoint()
    end
end

function _default_precondition(A, ::HansenBliekRohn)
    if is_H_matrix(A)
        return NoPrecondition()
    else
        return InverseMidpoint()
    end
end

# fallback
_default_precondition(_, ::AbstractLinearSolver) = InverseMidpoint()
