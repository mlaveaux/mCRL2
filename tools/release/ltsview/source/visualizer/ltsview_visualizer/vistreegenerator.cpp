#include "vistreegenerator.h"
#include "cluster.h"
#include "mcrl2/utilities/logger.h"
#include "glvistree.h"



void ConeConvertFunctor::matrixScale(float x, float y, float z){
    current_matrix.scale(x, y, z);
}

void ConeConvertFunctor::matrixTranslate(float x, float y, float z){
    current_matrix.translate(x, y, z);
}

void ConeConvertFunctor::matrixTranslate(const QVector3D& v){
    current_matrix.translate(v);
}

void ConeConvertFunctor::matrixRotate(float angle, float x, float y, float z){
    current_matrix.rotate(angle, x, y, z);
}


VisTree::VisTreeNode* ConeConvertFunctor::operator()(VisTree::VisTreeNode* parent, Cluster* cluster){
    if (parent == nullptr) current_matrix.setToIdentity(); else current_matrix = parent->data.matrix;

    // move down and to the side (if applicable)
    float branchtilt = cluster->isCentered() ? 0 : Settings::instance().branchTilt.value();
    float sideways = cluster->isCentered() ? 0 : cluster->getAncestor()->getBaseRadius() - cosf(branchtilt)*Settings::instance().clusterHeight.value();
    float down = cluster->hasDescendants() ? -Settings::instance().clusterHeight.value() : 0;
    // TODO: Check order
    matrixRotate(cluster->getPosition(), 0, 0, 1);
    matrixTranslate(sideways, 0, 0);
    matrixRotate(-branchtilt, 0, 1, 0);
    matrixTranslate(0, 0, down);

    if (!cluster->hasDescendants()){
        VisTree::VisTreeNode* node = new VisTree::VisTreeNode();
        GlUtil::Shapes::Sphere* sphere = new GlUtil::Shapes::Sphere();
        sphere->radius = 1;
        // find shape matrix
        float r = cluster->getTopRadius();
        matrixScale(r, r, r);

        node->data.matrix = current_matrix;
        node->data.shape = sphere;
        
        std::cout << "Sphere at: " << qvec4ToString(current_matrix.column(3)) << " with:" <<
        "\n\t- radius: 1\n\t- r: " << std::to_string(r) << std::endl;
        return node;
    }

   

    VisTree::VisTreeNode* node = new VisTree::VisTreeNode();
    GlUtil::Shapes::TruncatedCone* cone = new GlUtil::Shapes::TruncatedCone();
    cone->height = Settings::instance().clusterHeight.value();
    /// TODO: Figure out why top-bot are sometimes inverted in the clustering algorithm
    // cone->fill_bot = !cluster->isCentered() || parent == nullptr;
    // cone->fill_top = cluster->getNumDescendants() > 1;
    cone->fill_top = !cluster->isCentered() || parent == nullptr;
    cone->fill_bot = cluster->getNumDescendants() > 1;
    node->data.matrix = current_matrix;
    cone->radius_top = cluster->getTopRadius();
    cone->radius_bot = cluster->getBaseRadius();
    node->data.shape = cone;
    std::cout << "Cone at: " << qvec4ToString(current_matrix.column(3)) << " with:" << 
    "\n\t- topclosed: " << std::to_string(cone->fill_top) << 
    "\n\t- botclosed: " << std::to_string(cone->fill_bot) <<
    "\n\t- topradius: " << std::to_string(cone->radius_top) <<
    "\n\t- botradius: " << std::to_string(cone->radius_bot) <<
    "\n\t- height:    " << std::to_string(cone->height) << std::endl;

    return node;
}


VisTree::VisTreeNode* TubeConvertFunctor::operator()(VisTree::VisTreeNode* parent, Cluster* cluster){
    std::cout << "TubeConvertFunctor: not implemented yet..." << std::endl;
    return nullptr;
}

VisTree::VisTreeNode* VisTreeGenerator::generateTubes(Cluster* root){
    return generateClusterTree<TubeConvertFunctor>(root);
}

VisTree::VisTreeNode* VisTreeGenerator::generateCones(Cluster* root){
    return generateClusterTree<ConeConvertFunctor>(root);
}

VisTree::VisTreeNode* VisTreeGenerator::generate(Mode mode, Cluster* root){
    switch (mode){
        case Mode::TUBES:
            return generateTubes(root);
            break;
        case Mode::CONES:
            return generateCones(root);
            break;
        default:
            mCRL2log(mcrl2::log::warning) << "Invalid mode selected in ltsview scenegraphgenerator. Defaulting to tubes." << std::endl;
            return generateCones(root);
    }
}




