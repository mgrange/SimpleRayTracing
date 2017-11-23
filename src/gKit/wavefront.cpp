
#include <ctype.h>
#include <cstdio>

#include "wavefront.h"


Mesh read_mesh( const char *filename )
{
    Mesh data= create_mesh(GL_TRIANGLES);
    
    FILE *in= fopen(filename, "rt");
    if(in == NULL)
    {
        printf("loading mesh '%s'... failed.\n", filename);
        return data;
    }
    
    printf("loading mesh '%s'...\n", filename);
    
    std::vector<vec3> positions;
    std::vector<vec2> texcoords;
    std::vector<vec3> normals;
    
    std::vector<int> idp;
    std::vector<int> idt;
    std::vector<int> idn;
    
    char line_buffer[1024];
    bool error= true;
    for(;;)
    {
        // charge une ligne du fichier
        if(fgets(line_buffer, sizeof(line_buffer), in) == NULL)
        {
            error= false;       // fin du fichier, pas d'erreur detectee
            break;
        }
        
        // force la fin de la ligne, au cas ou
        line_buffer[sizeof(line_buffer) -1]= 0;
        
        // saute les espaces en debut de ligne
        char *line= line_buffer;
        while(*line && isspace(*line))
            line++;
        
        if(line[0] == 'v')
        {
            float x, y, z;
            if(line[1] == ' ')          // position x y z
            {
                if(sscanf(line, "v %f %f %f", &x, &y, &z) != 3)
                    break;
                positions.push_back( make_vec3(x, y, z) );
            }
            else if(line[1] == 'n')     // normal x y z
            {
                if(sscanf(line, "vn %f %f %f", &x, &y, &z) != 3)
                    break;
                normals.push_back( make_vec3(x, y, z) );
            }
            else if(line[1] == 't')     // texcoord x y
            {
                if(sscanf(line, "vt %f %f", &x, &y) != 2)
                    break;
                texcoords.push_back( make_vec2(x, y) );
            }
        }
        
        else if(line[0] == 'f')         // triangle a b c, les sommets sont numerotes a partir de 1 ou de la fin du tableau (< 0)
        {
            idp.clear();
            idt.clear();
            idn.clear();
            
            int next;
            for(line= line +1; ; line= line + next)
            {
                idp.push_back(0); 
                idt.push_back(0); 
                idn.push_back(0);         // 0: invalid index
                
                next= 0;
                if(sscanf(line, " %d/%d/%d %n", &idp.back(), &idt.back(), &idn.back(), &next) == 3) 
                    continue;
                else if(sscanf(line, " %d/%d %n", &idp.back(), &idt.back(), &next) == 2)
                    continue;
                else if(sscanf(line, " %d//%d %n", &idp.back(), &idn.back(), &next) == 2)
                    continue;
                else if(sscanf(line, " %d %n", &idp.back(), &next) == 1)
                    continue;
                else if(next == 0)      // fin de ligne
                    break;
            }
            
            for(int v= 2; v +1 < (int) idp.size(); v++)
            {
                int idv[3]= { 0, v -1, v };
                for(int i= 0; i < 3; i++)
                {
                    int k= idv[i];
                    int p= (idp[k] < 0) ? (int) positions.size() + idp[k] : idp[k] -1;
                    int t= (idt[k] < 0) ? (int) texcoords.size() + idt[k] : idt[k] -1;
                    int n= (idn[k] < 0) ? (int) normals.size()   + idn[k] : idn[k] -1;
                    
                    if(t >= 0) vertex_texcoord(data, texcoords[t]);
                    if(n >= 0) vertex_normal(data, normals[n]);
                    
                    if(p < 0) break; // error
                    push_vertex(data, positions[p]);
                }
            }
        }        
    }
    
    fclose(in);
    
    if(error)
        printf("loading mesh '%s'...\n[error]\n%s\n\n", filename, line_buffer);
    
    return data;
}

int write_mesh( const Mesh& mesh, const char *filename )
{
    if(mesh.primitives != GL_TRIANGLES)
        return -1;
    if(mesh.positions.size() == 0)
        return -1;
    if(filename == NULL)
        return -1;
    
    FILE *out= fopen(filename, "wt");
    if(out == NULL)
        return -1;
    
    printf("writing mesh '%s'...\n", filename);
    
    for(unsigned int i= 0; i < (unsigned int) mesh.positions.size(); i++)
        fprintf(out, "v %f %f %f\n", mesh.positions[i].x, mesh.positions[i].y, mesh.positions[i].z);
    fprintf(out, "\n");
    
    bool has_texcoords= (mesh.texcoords.size() == mesh.positions.size());
    for(unsigned int i= 0; i < (unsigned int) mesh.texcoords.size(); i++)
        fprintf(out, "vt %f %f\n", mesh.texcoords[i].x, mesh.texcoords[i].y);
    fprintf(out, "\n");
    
    bool has_normals= (mesh.normals.size() == mesh.positions.size());
    for(unsigned int i= 0; i < (unsigned int) mesh.normals.size(); i++)
        fprintf(out, "vn %f %f %f\n", mesh.normals[i].x, mesh.normals[i].y, mesh.normals[i].z);
    fprintf(out, "\n");
    
    bool has_indices= (mesh.indices.size() > 0);
    unsigned int n= has_indices ? (unsigned int) mesh.indices.size() : (unsigned int) mesh.positions.size();
    
    for(unsigned int i= 0; i +2 < n; i+= 3)
    {
        fprintf(out, "f");
        for(unsigned int k= 0; k < 3; k++)
        {
            unsigned int id= has_indices ? mesh.indices[i+k] +1 : i+k +1;
            fprintf(out, " %u", id);
            if(has_texcoords && has_normals)
                fprintf(out, "/%u/%u", id, id);
            else if(has_texcoords)
                fprintf(out, "/%u", id);
            else if(has_normals)
                fprintf(out, "//%u", id);
        }
        fprintf(out, "\n");
    }
    
    fclose(out);
    return 0;
}