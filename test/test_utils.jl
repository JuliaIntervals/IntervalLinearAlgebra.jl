@testset "Utils" begin
    orthants = Orthants(3)

    @test length(orthants) == 8
    eltype(Orthants) == Vector{Int}
    @test first(orthants) == [1, 1, 1]
    @test last(orthants) == [-1, -1, -1]
    @test orthants[4] == [-1, -1, 1]
end
