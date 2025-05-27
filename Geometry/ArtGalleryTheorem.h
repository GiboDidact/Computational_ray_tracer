#pragma once
#include "Polygon.h"

/*
Theorem 1.32: To cover a polygon with n vertices, floor(n/3) guards are needed for some polygons, and sufficient for all of them

This is a proof by Steve Fisk which uses Node coloring theorem. Proof is pretty simple, if you can just triangulate the polygon and for 
each triangle give the vertices 3 different colors then you just place a guard in one of the colors, since every triangle has 1 of each color
it can see the whole polygon. Its just showing this 3-color partition is possible

its just induction for a triangle n=3 its trivially true. Now if we have a polygon size n well an ear exists, so take the vertex at the tip of the ear out
and you get a polygon n-1 size which by I.H. is 3-colorable. Now add the tip back and just make it the color that isn't its 2 neighbors. Therefore the polygon is
3-colorable. It partitions it into 3 groups evenly distributed so get the minimum one floor(N/3) and you are done



Finding the Visibility Polygon for g. V(g) is the set of all points in the polygon that are visible to g. This is a polygon. How do we find this???

for visibility there seems to be "caves", theres a reflex vertex that makes it not visible and eventually the other opening will become visible.
no reflex vertex no hidden walls. At your point g you just spin around and travel along the boundary in a montonically increasing angle. If the polygon takes an opposite direction ignore it
and just find the intersection so you can keep going the same direction, do that until you reach the start


Winding Theorem: if you have a simply closed polygon, if you take a point inside and travel along all vertices it will form 1 CW or CCW loop. 
Has to do with orientation and being inside the shape or outside. Not sure of the proof but this has to be true.

Cave Theorem: If you are at a point G in a polygon and you have a visible vertex R, you have the segment GR, if walking to the next vertex is an "outside" turn,
if instead you kept walking forward (extending segment GR) you would eventually hit the boundary of the polygon.
(depending on whether its counterclockwise or clockwise oriented pretty much if you are going clockwise a left turn is outside, CCW a right turn is outside)
*intuitively its just saying whenever you have these outside reflex turns in a polygon it always has to "come back" or else it wouldnt be closed,
pretty much a consequence of the winding theorem, almost the same thing.

proof:the point G and R are inside the polygon and the segment is a diangol so its all inside the polygon too. By jordan curve theorem if a path doesnt intersect
the boundary than all points on a path are either completey inside or outside. Firstly, we know the segment GR doesnt end at vertex R. This is because
the Winding Theorem, since we are traveling CCW, WLOG, around the polygon it must form 1 CW loop at completion. Since we take a right turn we will have 
to take a left turn back to whatever angle our segment GR was at, or else how would it ever form a CCW loop? So we can extend the segment and it will hit another part of the
boundary by jordan theorem (if we didnt hit a boundary that means polygon is infinite and every point in the infinite line is inside the polygon).


Generate Visiblity Polygon V(G) Algorithm:

1.) given point, G, we know its in a triangle, pick a vertex in that triangle. We know that vertex is visible to G, add that to V(G)
2.) now go in order to the next vertex
  a.) if its visible then add that to V(G). if both vertices are visible to G then the whole edge is visible
  b.) if its not visible 2 cases arise (CCW outside turn is right, CW outside turn is left)
    c.) you took an outside turn (G,A-1,A): since you took an outside turn by the Cave Theorem we know that we can extend the segment G(A-1) and find a closest 
    intersection with the boundary at the angle we left at. Get this intersection point and add to V(G). We know all the previous vertices we skipped are 
    not visible to G because the "outside" step blocks all the "outside" vertices, and our intersection is the first one that becomes visible again.
    Since we took a right turn and then a left again you can imagine the boundaries crossing over the same line of sight angle, nothing would be visible
    to us. Find the next vertex of P on the edge you hit and continue.

    d.) you took an inside turn(G,A-1,A): the segment GA is blocked not because an "outside" turn but because another vertex/edge is blocking it. This is pretty
    much the opposite case of c now where at that "intersection point" and need to find the "outside" turn vertex. Pretty much you just walk along each vertex
    until you find the first visible vertex B. Thats the one blocking it and has to exist or else it wouldn't be blocked. Do the intersection again, 
    add to V(G) and continue from B. This vertex has to exist because G(A-1) is visible and as you travel along (A-1)A it gets blocked.
    The first thing to block it has to be a vertex or else an edge wouldve already blocked it.
    All the previous vertices skipped are not visible because we stopped at the first visible one.

3.) once we've gotten back to the start we are done. We winded around once and along every angle we found a place on the boundary. Therefore the 
Visibility Polygon, V(G), is correct because at every angle from G it contains the segment from G all the way to the visible boundary. And all the vertices
added to V(G) are inorder from P, therefore V(G) itself is in order and forms a closed curve and hence a polygon. 

Not that optimal in a way, it kinda is just ray casting, it uses a lot of intersections to see if points are visible, but is a little smarter about it,
its pretty much N^2



The smarter algorithms I think keep track of the angles and can tell if its visible as of consequence:
https://cs.uwaterloo.ca/research/tr/1985/CS-85-38.pdf

https://arxiv.org/pdf/1403.3905
*/

namespace ArtGalleryTheorem
{
    static void FindIntersectionPointOnBoundary(Polygon& P, glm::vec2 pos, int vertex_id, glm::ivec2& intersection_vertex, glm::vec2& intersection_pos, 
        bool know_edge = false, glm::ivec2 original_edge = glm::ivec2(0,0))
    {
        //create a ray (pos, P.V()[vertex_id] - pos) and get the first intersection on the boundary
        glm::vec2 segment_extend = pos + (10000.0f * glm::normalize(P.V()[vertex_id] - pos));
        float closest_distance = 1000000000.0f;
        for (int i = 0; i < P.V().size(); i++)
        {
            int e1 = i;
            int e2 = (i + 1) % P.V().size();
            if (know_edge && (e1 != original_edge.x || e2 != original_edge.y))
            {
                continue;
            }
            if (e1 != vertex_id && e2 != vertex_id)
            {
                //handle close to endpoint precision issues probably
                glm::vec3 result = Segment::LineToLineIntersectionPoint2(pos, segment_extend, P.V()[e1], P.V()[e2]);
                if (result.x == 1)
                {
                    float distance = glm::distance(pos, glm::vec2(result.y, result.z));
                    //std::cout << "distance: " << distance << std::endl;
                    if (distance < closest_distance)
                    {
                        closest_distance = distance;
                        intersection_vertex = glm::ivec2(e1, e2);
                        intersection_pos = glm::vec2(result.y, result.z);
                    }
                }
                else
                {
                    if (result.y == 999998)
                        std::cout << "lines parallel\n";
                }
            }
        }
        if (closest_distance == 1000000000.0f)
        {
            std::cout << "couldn't find an intersection??\n";
            std::cout << glm::normalize(P.V()[vertex_id] - pos).x<<" "<< glm::normalize(P.V()[vertex_id] - pos).y;
        }
    }

    static bool IsSegmentVisible(Polygon& P, glm::vec2 pos, int next_id)
    {
        //see if our position is right on top of the vertex segment would be like 0
        float eps2 = 0.0001f;
        if (P.hitVertex(pos, eps2) == 1)//next_id???
            return true;

        //see if its visible or not
        bool intersection = false;
        for (int i = 0; i < P.V().size(); i++)
        {
            int e1 = i;
            int e2 = (i + 1) % P.V().size();

            if (Segment::SegmentsIntersectNoBoundary(pos, P.V()[next_id], P.V()[e1], P.V()[e2]))
            {
                intersection = true;
                break;
            }
        }
        return !intersection;
    }


    //pretty much I brute force the visibility check which is n^2. A smarter method would use some wind sweeping method that keeps track of angles and such
    static Polygon GenerateVisibilityPolygon(Polygon P, glm::vec2 pos)
    {
        Timer time;
        time.Begin();
        // std::cout << "----GenerateVisibilityPolygon----\n";
        Polygon poly;

        if (!P.IsInside(pos))
            return poly;

        if (P.IsConvex())
            return P;
      
        //calculate visibility of every vertex from pos
        std::vector<bool> visible_vertex;
        int visible_id = -1;
        for (int i = 0; i < P.V().size(); i++)
        {
            bool visibility = IsSegmentVisible(P, pos, i);
            visible_vertex.push_back(visibility);
            if (visibility == true)
                visible_id = i;
        }
        if (visible_id == -1)
            std::cout << "no visible diagnols?\n";
       // std::cout << "time 1: " << time.getTimeMicro() << std::endl;
        //std::cout << "triangle id's: " << "(" << P.T()[triangle_id].x << ", " << P.T()[triangle_id].y << ", " << P.T()[triangle_id].z << ")\n";


        int vertex_id = visible_id;
        int starting_id = vertex_id;
        poly.AddVertex(P.V()[vertex_id]);
        //std::cout << "starting vertex: " << vertex_id << std::endl;

        int outside_turn_val = (P.IsCCW()) ? 1 : -1;
        int vertex_steps = 0;
        while (vertex_id != starting_id || vertex_steps == 0)
        {
            int next_id = (vertex_id + 1) % P.V().size();           
            bool visible = visible_vertex[next_id];

            if (visible)
            {
                poly.AddVertex(P.V()[next_id]);
                //std::cout << "visible\n";
            }
            else
            {
                //std::cout << "not visible\n";
                int turn_val = Segment::TurnSegments(pos, P.V()[vertex_id], P.V()[next_id]);
               // std::cout << turn_val << std::endl;
                if (turn_val == outside_turn_val)
                {

                    //turned outside
                    //std::cout << "outside turn\n";
                    glm::ivec2 intersection_vertex;
                    glm::vec2  intersection_pos;
                    FindIntersectionPointOnBoundary(P, pos, vertex_id, intersection_vertex, intersection_pos);

                    poly.AddVertex(intersection_pos);
                    //set it to the one before so next loop we step to the one we are supposed to
                    vertex_id = intersection_vertex.x;
                    
                    vertex_steps++;
                    continue;
                }
                else if(turn_val != 0)
                {
                    //turned inside
                    //std::cout << "inside turn\n";
                    //go until you find the first visible vertex
                    bool visible = false;  
                    glm::ivec2 original_edge(vertex_id, next_id);
                    while (!visible)
                    {
                        next_id = (next_id + 1) % P.V().size();
                        visible = visible_vertex[next_id];
                    }

                    //next_id is our vertex we want to do query with, do intersection
                    glm::ivec2 intersection_vertex;
                    glm::vec2  intersection_pos;
                    //*we know what edge we have to intersect 
                    FindIntersectionPointOnBoundary(P, pos, next_id, intersection_vertex, intersection_pos, true, original_edge);
                    
                    //set vertex_id to the visible vertex we hit
                    poly.AddVertex(intersection_pos);
                    if(next_id != starting_id)
                        poly.AddVertex(P.V()[next_id]);
                    vertex_id = next_id;

                    vertex_steps++;
                    continue;
                }
                else if (turn_val == 0)
                {
                    //its colinear, i guess its visible
                    poly.AddVertex(P.V()[next_id]);
                }
            }

            vertex_id = next_id;
            vertex_steps++;
        }
       // std::cout << "finish time: " << time.getTimeMicro() << std::endl;

        if (!poly.isValidPolygon())
            std::cout << "not a valid polygon!\n";
        poly.TriangulateDiagnolSplitting(true);
        //std::cout << "vertex steps: " << vertex_steps <<"polygon size: "<<poly.V().size()<< std::endl;
        return poly;
    }

    static void PlaceGuardsColorTheoremFindEar(Polygon P, std::vector<glm::ivec3>& ears, Polygon& original_P)
    {
        //find the ear, then mark the triangle as found and take note of ear
        //then for earintriangulation just pass in marked triangles and ignore
        if (P.V().size() == 3)
        {
            glm::ivec3 converted_ids = glm::ivec3(original_P.findVertexIndex(P.V()[0]), original_P.findVertexIndex(P.V()[1]),
                original_P.findVertexIndex(P.V()[2]));

            ears.push_back(converted_ids);

            return;
        }

        glm::ivec4 results = P.FindEarInTriangulation();
        glm::ivec3 ear(results.x, results.y, results.z);
        if (ear.x == -1 && ear.y == -1) {
            std::cout << "no ears?\n";
            return;
        }

        glm::ivec3 converted_ids = glm::ivec3(original_P.findVertexIndex(P.V()[ear.x]), original_P.findVertexIndex(P.V()[ear.y]),
                                              original_P.findVertexIndex(P.V()[ear.z]));
        
        ears.push_back(converted_ids);
        P.RemoveVertexEar(ear.y);
        PlaceGuardsColorTheoremFindEar(P, ears, original_P);
    }

    //worst case: floor(n/3) guards
    //so just find an ear, place it on the array in order, and take the tip of the ear out of the polygon and repeat again until the polygon is
    //size 3. Then all you do is give the triangle 3 colors, and reconstruct one point by one point adding ear tips and picking the other color
    static std::vector<glm::vec2> PlaceGuardsColorTheorem(Polygon P)
    {
        //find ears
        std::vector<glm::ivec3> ears;
        PlaceGuardsColorTheoremFindEar(P, ears, P);
        std::reverse(ears.begin(), ears.end());
        
        //give colors to the triangle ,-1 is default, 0, 1 ,2
        std::vector<int> colors(P.V().size(), -1); //maps P.V()-->color
        colors[ears[0].x] = 0;
        colors[ears[0].y] = 1;
        colors[ears[0].z] = 2;

        for (int i = 1; i < ears.size(); i++)
        {
            //eartip is the middle one, so just find a different color
            if (colors[ears[i].y] != -1)
                std::cout << "eartip already colored??\n";

            if((colors[ears[i].x] == 1 && colors[ears[i].z] == 2) || ((colors[ears[i].x] == 2 && colors[ears[i].z] == 1)))
                colors[ears[i].y] = 0;
            else if ((colors[ears[i].x] == 0 && colors[ears[i].z] == 2) || ((colors[ears[i].x] == 2 && colors[ears[i].z] == 0)))
                colors[ears[i].y] = 1;
            else if ((colors[ears[i].x] == 0 && colors[ears[i].z] == 1) || ((colors[ears[i].x] == 1 && colors[ears[i].z] == 0)))
                colors[ears[i].y] = 2;
        }

        
        //now place guards at every vertex that has color 0
        std::vector<glm::vec2>  guards;
        int zero_counter = 0;
        int one_counter = 0;
        int two_counter = 0;
        int lowest_color = -1;
        for (int i = 0; i < colors.size(); i++)
        {
            if (colors[i] == 0)
                zero_counter++;
            if (colors[i] == 1)
                one_counter++;
            if (colors[i] == 2)
                two_counter++;
            //std::cout << colors[i] << ", ";
        }
        if(zero_counter <= one_counter && zero_counter <=two_counter)
            lowest_color = 0;
        else if (one_counter <= zero_counter && one_counter <= two_counter)
            lowest_color = 1;
        else if (two_counter <= one_counter && two_counter <= zero_counter)
            lowest_color = 2;

        for (int i = 0; i < colors.size(); i++)
        {
            if (colors[i] == lowest_color)
            {
                glm::vec2 pos1 = P.V()[i];

                 //cant put it exactly on the vertex so nudge it inside a little
                 int previous = (i + colors.size() - 1) % colors.size();
                 int forward = (i + 1) % colors.size();
                 glm::vec2 dir = P.V()[forward] - P.V()[previous];
                 if(P.isVertexReflex(i))
                    dir = (glm::normalize((P.V()[i] - P.V()[forward])) + glm::normalize((P.V()[i] - P.V()[previous]))) / 2.0f;
                 else
                    dir = (glm::normalize((P.V()[forward]- P.V()[i])) + glm::normalize((P.V()[previous])- P.V()[i])) / 2.0f;
     
                 dir = glm::normalize(dir);
                 float eps_t = 0.01f; 
                 //have to figure out if t is positive or negative, just figure out if it a left or right turn
                 glm::vec2 nudge_1 = pos1 + dir * eps_t;
                 glm::vec2 nudge_2 = pos1 - dir * eps_t;
                 int outside_turn_val = (P.IsCCW()) ? 1 : -1;
                 if (Segment::TurnSegments(P.V()[previous], P.V()[i], nudge_1) != outside_turn_val)
                 {
                     guards.push_back(nudge_1);
                 }
                 else
                 {
                     guards.push_back(nudge_2);
                 }
            }
        }
       // std::cout << std::endl;
       // std::cout << floor(P.V().size()/3.0f) <<"  "<< zero_counter << " " << one_counter << " " << two_counter << std::endl;
        
        return guards;
    }

    //n-2 guards
    static std::vector<glm::vec2> PlaceGuardsTriangles(Polygon P)
    {
        std::vector<glm::vec2>  guards;
        for (int i = 0; i < P.T().size(); i++)
        {
            glm::vec2 pos1 = P.V()[P.T()[i].x];
            glm::vec2 pos2 = P.V()[P.T()[i].y];
            glm::vec2 pos3 = P.V()[P.T()[i].z];

            glm::vec2 center = (pos1 + pos2 + pos3) / 3.0f;

            guards.push_back(center);
        }
        return guards;
    }
    
    static int TriangleShareEdge(Polygon& P, int tri_ignore, int e1, int e2)
    {
        for (int t = 0; t < P.T().size(); t++)
        {
            if (t == tri_ignore)
                continue;
            if ((P.T()[t].x == e1 || P.T()[t].y == e1 || P.T()[t].z == e1) && (P.T()[t].x == e2 || P.T()[t].y == e2 || P.T()[t].z == e2))
            {
                return t;
            }
        }
        return -1;
    }

    //n/2 guards
    static std::vector<glm::vec2> PlaceGuardsTrianglesShared(Polygon P)
    {
        std::vector<glm::vec2>  guards;

        std::vector<bool> taken(P.T().size(), false);
        for (int i = 0; i < P.T().size(); i++)
        {
            //give your edge and find triangle that shares that edge
            int v1 = P.T()[i].x;
            int v2 = P.T()[i].y;
            int v3 = P.T()[i].z;

            int n1 = TriangleShareEdge(P, i, v1, v2);
            int n2 = TriangleShareEdge(P, i, v2, v3);
            int n3 = TriangleShareEdge(P, i, v1, v3);

            if (n1 != -1 && !taken[n1])
            {
                glm::vec2 pos1 = P.V()[v1];
                glm::vec2 pos2 = P.V()[v2];

                glm::vec2 center = (pos1 + pos2) / 2.0f;
                guards.push_back(center);
                taken[n1] = true;
                taken[i] = true;
            }
            else if (n2 != -1 && !taken[n2])
            {
                glm::vec2 pos1 = P.V()[v2];
                glm::vec2 pos2 = P.V()[v3];

                glm::vec2 center = (pos1 + pos2) / 2.0f;
                guards.push_back(center);
                taken[n2] = true;
                taken[i] = true;
            }
            else if (n3 != -1 && !taken[n3])
            {
                glm::vec2 pos1 = P.V()[v1];
                glm::vec2 pos2 = P.V()[v3];

                glm::vec2 center = (pos1 + pos2) / 2.0f;
                guards.push_back(center);
                taken[n3] = true;
                taken[i] = true;
            }
            else if(!taken[i])
            {
                glm::vec2 pos1 = P.V()[P.T()[i].x];
                glm::vec2 pos2 = P.V()[P.T()[i].y];
                glm::vec2 pos3 = P.V()[P.T()[i].z];
                
                glm::vec2 center = (pos1 + pos2 + pos3) / 3.0f;
                guards.push_back(center);
                taken[i] = true;
            }
        }
        return guards;
    }

    //another method is rank every vertex and prioritize the higher rank ones
    //or just do the one method and throw away guards that can see eachother
}