#pragma once

#include "gene/core/dimension.hpp"
#include "onyx/app/app.hpp"
#include "onyx/rendering/render_context.hpp"

namespace Gene::Demo
{
template <Dimension D> class Layer : public Onyx::UserLayer
{
  public:
    Layer(Onyx::SingleWindowApp *p_Application);

    void OnUpdate() override;

  private:
    Onyx::SingleWindowApp *m_Application;
    Onyx::Window *m_Window;
    Onyx::RenderContext<D> *m_Context;
    Onyx::Camera<D> *m_Camera;
};
} // namespace Gene::Demo
