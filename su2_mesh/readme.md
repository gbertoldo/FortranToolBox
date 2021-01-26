Este programa gera um script para geração de malhas estruturadas e não estruturadas com o programa Gmsh.
As malhas geradas servem para otimização de cones nasais em regime supersônico.

Para compilação, execute o seguinte comando:
./compile.sh

O script gera um executável denominado su2mesh.x.

Edite o arquivo input.txt para definir os parâmetros de configuração da malha.

Execute o programa:
./su2mesh.x

O programa cria um script denominado mesh.geo. Abra o script com o Gmsh para a geração da malha.
