function solve(A, b,
               solver::IterativeSolver,
               precondition::Precondition=_default_precondition(A, solver))

    A, b = precondition(A, b)
    x = enclose(A, b)

    solver(x, A, b)

    return x
end

function solve(A, b,
               solver::DirectSolver,
               precondition::Precondition=_default_precondition(A, solver))

    A, b = precondition(A, b)
    return solver(A, b)
end


# fallback
function solve(A, b,
               solver::LinearSolver=_default_solver(),
               precondition::Precondition=_default_precondition(A, solver))

    A, b = precondition(A, b)
    return solver(A, b)
end

## Default settings
_default_solver() = GaussElimination()

function _default_precondition(A, ::GaussElimination)
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
_default_precondition(_, ::LinearSolver) = InverseMidpoint()
