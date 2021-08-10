@testset "Utils" begin
    @test list_orthants(2) == [[1, 1], [-1, 1], [-1, -1], [1, -1]]
end
